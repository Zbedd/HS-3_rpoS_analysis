"""Per-sector colony detection.

Default detection: per-wedge median reference + 96th-percentile
distance threshold + small morphology + connected-component selection
+ post-CC smoothing close + flood-fill hole fill. Robust across the
full set of plates and naturally accommodates dim / dark colonies
(e.g. rpoS) that read as low-contrast against agar.

Three alternative methods are available for stubborn cases:
  - `blank_otsu`  : reference = blank-sector median, threshold = Otsu
                    on wedge distance histogram (with sanity cap)
  - `raw_p96`     : default logic on the raw linear-light image
                    (rescues dim R2A colonies that calibration crushes
                    into noise)
  - `convex_hull` : runs raw_p96 then takes the convex hull of the
                    contour (rescues colonies whose rim is the only
                    high-contrast feature)

`detect_colonies` runs the default everywhere, then applies any
manual `TARGETED_ADJUSTMENTS`, then auto-overrides any remaining
spots with low solidity / small area.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import cv2
import numpy as np

from .config import (
    BLANK_SECTOR, NA_SPOTS, N_SECTORS, SECTOR_0_THETA_DEG, SECTOR_SAMPLE,
    SPOT_AREA_RATIO_FLOOR, SPOT_CENTER_ANGLE_TOL_DEG, SPOT_DIST_PERCENTILE,
    SPOT_MAX_AREA_FRAC, SPOT_MIN_AREA_FRAC, SPOT_MIN_AREA_PX,
    SPOT_MORPH_KSIZE, SPOT_OTSU_FALLBACK_PERCENTILE, SPOT_OTSU_MAX_FRAC,
    SPOT_R_INNER_FRAC, SPOT_R_OUTER_FRAC, SPOT_SMOOTH_KSIZE,
    SPOT_SOLIDITY_FLOOR, TARGETED_ADJUSTMENTS, PLATE_INTERIOR_FRAC,
)
from .geometry import (
    Geometry, _angle_grid, _radius_grid, _sector_mask,
)


# ── Spot dataclass ─────────────────────────────────────────────

@dataclass
class Spot:
    sample: str
    sector: int
    cx: float
    cy: float
    r_eff: float
    area_px: int
    contour: list  # list of [x, y] ints (CV-style outer contour)
    below_floor: bool = False  # area_px < SPOT_MIN_AREA_PX
    solidity: float = 1.0      # contour_area / convex_hull_area
    auto_override: str = ''    # method used if auto-overridden, else ''


def _solidity(contour) -> float:
    """contour-fill area / convex-hull area. 1.0 == convex shape; ~0.5
    for a crescent. Returns 0.0 for degenerate contours."""
    pts = np.asarray(contour, dtype=np.int32)
    if pts.shape[0] < 3:
        return 0.0
    a = float(cv2.contourArea(pts))
    a_hull = float(cv2.contourArea(cv2.convexHull(pts)))
    if a_hull <= 0.0:
        return 0.0
    return a / a_hull


# ── Common helpers ─────────────────────────────────────────────

def _fill_holes(mask_u8: np.ndarray) -> np.ndarray:
    """Fill internal holes in a uint8 binary mask. Standard OpenCV
    idiom: flood-fill the background from the image corner and OR the
    inverse — any 0-pixel not reached by the flood is an internal
    hole.

    Handles true donuts (closed inner ring). C-shapes — rings with a
    thin outward gap — are NOT closed by this; those are addressed
    upstream by the blank-sector reference (which prevents the colony
    centre from reading as background)."""
    h, w = mask_u8.shape
    flood = mask_u8.copy()
    ff_mask = np.zeros((h + 2, w + 2), dtype=np.uint8)
    cv2.floodFill(flood, ff_mask, (0, 0), 255)
    holes = cv2.bitwise_not(flood)
    return cv2.bitwise_or(mask_u8, holes)


def _finalise_component(comp_mask: np.ndarray, sector: int) -> Optional[dict]:
    """Common post-CC pipeline: smooth (CLOSE only), fill internal
    holes, re-derive area / centroid, extract the outer contour."""
    k_smooth = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (SPOT_SMOOTH_KSIZE, SPOT_SMOOTH_KSIZE))
    comp_mask = cv2.morphologyEx(comp_mask, cv2.MORPH_CLOSE, k_smooth)
    comp_mask = _fill_holes(comp_mask)

    ys, xs = np.where(comp_mask > 0)
    if xs.size == 0:
        return None
    area = int(xs.size)
    cx = float(xs.mean())
    cy = float(ys.mean())

    contours, _ = cv2.findContours(comp_mask, cv2.RETR_EXTERNAL,
                                    cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        raise RuntimeError(f'failed to extract contour for sector {sector}')
    contour = max(contours, key=cv2.contourArea).reshape(-1, 2).astype(int).tolist()
    return {'sector': sector, 'cx': cx, 'cy': cy,
            'area': area, 'contour': contour,
            'below_floor': area < SPOT_MIN_AREA_PX}


def _select_cc(binary: np.ndarray, geom: Geometry, sector: int,
               min_area: float, max_area: float) -> Optional[np.ndarray]:
    """Pick the candidate component for `sector`: the largest CC
    within `SPOT_CENTER_ANGLE_TOL_DEG` of the wedge centerline that
    satisfies the area band. Returns the chosen CC mask (uint8) or
    None. Ranking by **area first, angle as tiebreaker** keeps a
    small dish-edge artefact at the centerline from beating a larger
    compact colony slightly off-axis."""
    span = 360.0 / N_SECTORS
    n_lab, lab, stats, centroids = cv2.connectedComponentsWithStats(
        binary, connectivity=8)
    wedge_center_theta = (SECTOR_0_THETA_DEG + sector * span) % 360.0
    candidates = []
    for i in range(1, n_lab):
        area = int(stats[i, cv2.CC_STAT_AREA])
        if area < min_area or area > max_area:
            continue
        cand_cx, cand_cy = centroids[i]
        cand_theta = float(np.degrees(np.arctan2(
            -(cand_cy - geom.dish_cy), cand_cx - geom.dish_cx)))
        if cand_theta < 0:
            cand_theta += 360.0
        angle_off = abs((cand_theta - wedge_center_theta + 180.0) % 360.0 - 180.0)
        if angle_off > SPOT_CENTER_ANGLE_TOL_DEG:
            continue
        candidates.append((-area, angle_off, i))
    if not candidates:
        return None
    candidates.sort()
    label_id = candidates[0][2]
    return (lab == label_id).astype(np.uint8) * 255


# ── Per-sector detectors ───────────────────────────────────────

def _detect_one_sector(img: np.ndarray, geom: Geometry, sector: int,
                       interior: np.ndarray, theta: np.ndarray,
                       min_area: float, max_area: float) -> Optional[dict]:
    """Default detection: per-wedge median reference + p96 distance
    threshold."""
    wedge = interior & _sector_mask(theta, sector)
    if not wedge.any():
        raise RuntimeError(f'sector {sector} has no in-band pixels')
    h, w = img.shape[:2]
    wedge_pix = img[wedge]
    med = np.median(wedge_pix, axis=0).astype(np.float32)
    dist_wedge = np.linalg.norm(wedge_pix - med, axis=1)
    thresh = float(np.percentile(dist_wedge, SPOT_DIST_PERCENTILE))
    binary = np.zeros((h, w), dtype=np.uint8)
    binary[wedge] = (dist_wedge >= thresh).astype(np.uint8) * 255
    k = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (SPOT_MORPH_KSIZE, SPOT_MORPH_KSIZE))
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, k)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, k)
    comp_mask = _select_cc(binary, geom, sector, min_area, max_area)
    if comp_mask is None:
        return None
    return _finalise_component(comp_mask, sector)


def _detect_one_sector_raw_p96(img_raw: np.ndarray, geom: Geometry,
                               sector: int, interior: np.ndarray,
                               theta: np.ndarray,
                               min_area: float,
                               max_area: float) -> Optional[dict]:
    """Same logic as `_detect_one_sector` but on the raw linear-light
    image. Per-angle calibration anchors black ring → 0, which on
    plates where agar reflects darker than the printed black pushes
    the entire dish to near-zero negative values — colonies and agar
    squash into a tight, very dim band where the threshold pass
    can't separate them. Detection on raw preserves the original
    colour contrasts and reliably segments dim colonies as clean
    disks.
    """
    wedge = interior & _sector_mask(theta, sector)
    if not wedge.any():
        raise RuntimeError(f'sector {sector} has no in-band pixels')
    h, w = img_raw.shape[:2]
    wedge_pix = img_raw[wedge]
    med = np.median(wedge_pix, axis=0).astype(np.float32)
    dist = np.linalg.norm(wedge_pix - med, axis=1)
    thresh = float(np.percentile(dist, SPOT_DIST_PERCENTILE))
    binary = np.zeros((h, w), dtype=np.uint8)
    binary[wedge] = (dist >= thresh).astype(np.uint8) * 255
    k = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (SPOT_MORPH_KSIZE, SPOT_MORPH_KSIZE))
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, k)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, k)
    comp_mask = _select_cc(binary, geom, sector, min_area, max_area)
    if comp_mask is None:
        return None
    return _finalise_component(comp_mask, sector)


def _convex_hull_result(r: dict) -> dict:
    """Replace `r`'s contour with its convex hull, recomputing area /
    centroid from the filled hull. Solidity becomes 1.0 by
    construction. For colonies whose rim is the only high-contrast
    feature against agar, the rim segments as a crescent / donut
    and the hull recovers the underlying filled disk."""
    pts = np.asarray(r['contour'], dtype=np.int32)
    if pts.shape[0] < 3:
        return r
    hull = cv2.convexHull(pts).reshape(-1, 2).astype(int)
    x0, y0 = hull.min(axis=0)
    x1, y1 = hull.max(axis=0)
    pad = 4
    H = int(y1 - y0 + 2 * pad)
    W = int(x1 - x0 + 2 * pad)
    hull_local = hull - [x0 - pad, y0 - pad]
    mask = np.zeros((H, W), dtype=np.uint8)
    cv2.fillPoly(mask, [hull_local], 255)
    ys, xs = np.where(mask > 0)
    if xs.size == 0:
        return r
    area = int(xs.size)
    cx = float(xs.mean()) + (x0 - pad)
    cy = float(ys.mean()) + (y0 - pad)
    return {**r, 'cx': cx, 'cy': cy, 'area': area,
            'contour': hull.tolist()}


def _detect_one_sector_convex_hull(img_raw: np.ndarray, geom: Geometry,
                                   sector: int, interior: np.ndarray,
                                   theta: np.ndarray,
                                   min_area: float,
                                   max_area: float) -> Optional[dict]:
    r = _detect_one_sector_raw_p96(img_raw, geom, sector, interior,
                                   theta, min_area, max_area)
    if r is None:
        return None
    return _convex_hull_result(r)


def _detect_one_sector_blank_otsu(img: np.ndarray, geom: Geometry,
                                  sector: int, interior: np.ndarray,
                                  theta: np.ndarray, agar_ref: np.ndarray,
                                  min_area: float,
                                  max_area: float) -> Optional[dict]:
    """Distance reference is the blank-sector median (guaranteed
    colony-free). Threshold is Otsu's method on the wedge distance
    histogram, with a sanity cap that falls back to a percentile
    when Otsu over-selects (weak-contrast wedge with no clean
    bimodal split).
    """
    wedge = interior & _sector_mask(theta, sector)
    if not wedge.any():
        raise RuntimeError(f'sector {sector} has no in-band pixels')
    h, w = img.shape[:2]
    wedge_pix = img[wedge]
    dist_wedge = np.linalg.norm(wedge_pix - agar_ref, axis=1)
    dmax = float(dist_wedge.max())
    if dmax <= 0.0:
        return None
    dist_u8 = np.clip(dist_wedge / dmax * 255.0, 0, 255).astype(np.uint8)
    otsu_t, _ = cv2.threshold(dist_u8, 0, 255,
                              cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    thresh = (otsu_t / 255.0) * dmax
    if float((dist_wedge >= thresh).mean()) > SPOT_OTSU_MAX_FRAC:
        thresh = float(np.percentile(dist_wedge,
                                     SPOT_OTSU_FALLBACK_PERCENTILE))
    binary = np.zeros((h, w), dtype=np.uint8)
    binary[wedge] = (dist_wedge >= thresh).astype(np.uint8) * 255
    k = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (SPOT_MORPH_KSIZE, SPOT_MORPH_KSIZE))
    binary = cv2.morphologyEx(binary, cv2.MORPH_CLOSE, k)
    binary = cv2.morphologyEx(binary, cv2.MORPH_OPEN, k)
    comp_mask = _select_cc(binary, geom, sector, min_area, max_area)
    if comp_mask is None:
        return None
    return _finalise_component(comp_mask, sector)


# ── Plate-level driver + overrides ─────────────────────────────

def detect_colonies(img: np.ndarray, geom: Geometry,
                    stem: Optional[str] = None,
                    img_raw: Optional[np.ndarray] = None) -> list:
    """Detect colonies in the five non-blank sectors. Returns the list
    of `Spot`s.

    Default detection runs on the calibrated image. If `stem` matches
    an entry in `TARGETED_ADJUSTMENTS`, the listed samples are
    re-detected with the configured alternative method. Then the
    auto-override pass tries `blank_otsu` and (if `img_raw` is
    provided) `raw_p96` on any remaining spots with low solidity or
    suspiciously small area, keeping whichever variant has the
    highest solidity.

    Throws if a non-blank sector (not configured as NA) has no
    qualifying colony candidate.
    """
    h, w = img.shape[:2]
    theta = _angle_grid(geom.dish_cx, geom.dish_cy, h, w)
    rgrid = _radius_grid(geom.dish_cx, geom.dish_cy, h, w)
    inner = SPOT_R_INNER_FRAC * geom.r_dish
    outer = SPOT_R_OUTER_FRAC * geom.r_dish
    interior = (rgrid >= inner) & (rgrid <= outer) \
               & (rgrid < geom.r_dish * PLATE_INTERIOR_FRAC)
    plate_area = float(np.pi * (geom.r_dish ** 2))
    min_area = SPOT_MIN_AREA_FRAC * plate_area
    max_area = SPOT_MAX_AREA_FRAC * plate_area

    na_samples = NA_SPOTS.get(stem, set()) if stem else set()
    if na_samples:
        for s in na_samples:
            print(f'  [NA] {stem}/{s}: skipping detection (configured)',
                  flush=True)

    spots = []
    missing = []
    for sector in range(N_SECTORS):
        if sector == BLANK_SECTOR:
            continue
        sample = SECTOR_SAMPLE[sector]
        if sample in na_samples:
            continue
        r = _detect_one_sector(img, geom, sector, interior, theta,
                               min_area, max_area)
        if r is None:
            missing.append((sector, sample))
            continue
        spots.append(Spot(
            sample=sample, sector=sector,
            cx=r['cx'], cy=r['cy'],
            r_eff=float(np.sqrt(r['area'] / np.pi)),
            area_px=r['area'], contour=r['contour'],
            below_floor=bool(r.get('below_floor', False)),
            solidity=_solidity(r['contour']),
        ))
    if missing:
        raise RuntimeError(
            f'no colony candidate in non-blank sector(s) {missing}; '
            f'plate may be misorientated (blank should be at top, sector 0) '
            f'or a colony failed to grow.'
        )

    overrides = TARGETED_ADJUSTMENTS.get(stem, {}) if stem else {}
    if overrides:
        spots = _apply_targeted_adjustments(
            spots, overrides, img, img_raw, geom, interior, theta,
            min_area, max_area, stem)

    spots = _apply_auto_overrides(spots, overrides, img, img_raw,
                                  geom, interior, theta,
                                  min_area, max_area, stem)
    return spots


def _apply_targeted_adjustments(spots: list, overrides: dict,
                                img: np.ndarray,
                                img_raw: Optional[np.ndarray],
                                geom: Geometry,
                                interior: np.ndarray, theta: np.ndarray,
                                min_area: float, max_area: float,
                                stem: str) -> list:
    """Re-detect specific (sample) entries on a plate using the
    method configured in `TARGETED_ADJUSTMENTS[stem]`. Replaces the
    matching `Spot` in `spots`; non-overridden spots are unchanged."""
    blank_wedge = interior & _sector_mask(theta, BLANK_SECTOR)
    agar_ref = None
    out = list(spots)
    for sample, cfg in overrides.items():
        method = cfg.get('method')
        try:
            sector = next(s for s, name in SECTOR_SAMPLE.items() if name == sample)
        except StopIteration:
            raise RuntimeError(
                f'targeted_adjustment for {stem!r}: unknown sample {sample!r}')
        if method == 'blank_otsu':
            if agar_ref is None:
                if not blank_wedge.any():
                    raise RuntimeError(
                        f'targeted_adjustment ({stem}/{sample}): blank sector '
                        f'has no pixels in the search annulus')
                agar_ref = np.median(img[blank_wedge], axis=0).astype(np.float32)
            r = _detect_one_sector_blank_otsu(
                img, geom, sector, interior, theta, agar_ref,
                min_area, max_area)
        elif method == 'raw_p96':
            if img_raw is None:
                raise RuntimeError(
                    f'targeted_adjustment ({stem}/{sample}): method '
                    f'raw_p96 requires the raw image (img_raw=None)')
            r = _detect_one_sector_raw_p96(
                img_raw, geom, sector, interior, theta,
                min_area, max_area)
        elif method == 'convex_hull':
            if img_raw is None:
                raise RuntimeError(
                    f'targeted_adjustment ({stem}/{sample}): method '
                    f'convex_hull requires the raw image (img_raw=None)')
            r = _detect_one_sector_convex_hull(
                img_raw, geom, sector, interior, theta,
                min_area, max_area)
        else:
            raise RuntimeError(
                f'targeted_adjustment for {stem!r}: unknown method {method!r}')
        if r is None:
            raise RuntimeError(
                f'targeted_adjustment {stem}/{sample} (method={method}) '
                f'produced no candidate')
        new_spot = Spot(
            sample=sample, sector=sector,
            cx=r['cx'], cy=r['cy'],
            r_eff=float(np.sqrt(r['area'] / np.pi)),
            area_px=r['area'], contour=r['contour'],
            below_floor=bool(r.get('below_floor', False)),
            solidity=_solidity(r['contour']),
            auto_override=f'manual:{method}',
        )
        out = [new_spot if s.sample == sample else s for s in out]
        print(f'  [adj] {stem}/{sample}: method={method} '
              f'r_eff={new_spot.r_eff:.1f} area={new_spot.area_px}',
              flush=True)
    return out


def _apply_auto_overrides(spots, manual_overrides, img, img_raw, geom,
                          interior, theta, min_area, max_area, stem):
    """Cascade for stubborn detections: try `blank_otsu` (calibrated)
    and `raw_p96` (uncalibrated), keep the candidate with the highest
    solidity. Skips spots already in `manual_overrides`."""
    if not spots:
        return spots
    plate_median = float(np.median([s.area_px for s in spots]))
    blank_wedge = interior & _sector_mask(theta, BLANK_SECTOR)
    agar_ref_cache = [None]

    def _agar_ref():
        if agar_ref_cache[0] is None:
            if not blank_wedge.any():
                raise RuntimeError(
                    f'auto-override ({stem}): blank sector has no '
                    f'pixels in the search annulus')
            agar_ref_cache[0] = np.median(img[blank_wedge], axis=0).astype(np.float32)
        return agar_ref_cache[0]

    out = []
    for s in spots:
        if s.sample in manual_overrides:
            out.append(s)
            continue
        bad_solidity = s.solidity < SPOT_SOLIDITY_FLOOR
        bad_area = s.area_px < SPOT_AREA_RATIO_FLOOR * plate_median
        if not (bad_solidity or bad_area):
            out.append(s)
            continue

        cands = [('default', None, s.solidity, s.area_px,
                  list(s.contour), s.cx, s.cy)]
        try:
            r_bo = _detect_one_sector_blank_otsu(
                img, geom, s.sector, interior, theta, _agar_ref(),
                min_area, max_area)
        except RuntimeError:
            r_bo = None
        if r_bo is not None:
            cands.append(('blank_otsu', r_bo, _solidity(r_bo['contour']),
                          r_bo['area'], r_bo['contour'],
                          r_bo['cx'], r_bo['cy']))

        if img_raw is not None:
            try:
                r_raw = _detect_one_sector_raw_p96(
                    img_raw, geom, s.sector, interior, theta,
                    min_area, max_area)
            except RuntimeError:
                r_raw = None
            if r_raw is not None:
                cands.append(('raw_p96', r_raw,
                              _solidity(r_raw['contour']),
                              r_raw['area'], r_raw['contour'],
                              r_raw['cx'], r_raw['cy']))

        # Pick the candidate with the highest solidity. Tie-break
        # by area in the direction that fixes the trigger.
        def score(c):
            label, _r, sol, area, _ct, _cx, _cy = c
            return (sol,
                    area if bad_area else -abs(area - plate_median))
        best = max(cands, key=score)
        label, r_chosen, sol_chosen, area_chosen, ct_chosen, cx_chosen, cy_chosen = best
        if label == 'default':
            out.append(s)
            continue

        # Only keep the override if it materially improves solidity
        # (or rescues a tiny rim into a much larger disk at
        # comparable solidity). Avoids swapping in a noise-driven
        # retry.
        keep = (sol_chosen > s.solidity + 0.05) or (
            bad_area and area_chosen >= 1.5 * s.area_px
            and sol_chosen >= s.solidity)
        if not keep:
            out.append(s)
            continue

        out.append(Spot(
            sample=s.sample, sector=s.sector,
            cx=cx_chosen, cy=cy_chosen,
            r_eff=float(np.sqrt(area_chosen / np.pi)),
            area_px=area_chosen, contour=ct_chosen,
            below_floor=bool(r_chosen.get('below_floor', False)),
            solidity=sol_chosen,
            auto_override=label,
        ))
        print(f'  [auto] {stem}/{s.sample}: method={label} '
              f'solidity {s.solidity:.2f}->{sol_chosen:.2f}, '
              f'area {s.area_px}->{area_chosen} '
              f'(reason: {"crescent" if bad_solidity else "small"})',
              flush=True)
    return out
