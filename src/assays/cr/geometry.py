"""Plate geometry: red-ring centre + radii, dish boundary, sectoring.

The geometry pipeline runs in two stages:
  1. `_detect_red_ring` finds the printed red ring's centre and a
     spec-derived inner / outer radius from a Hough fit on the
     strict-red mask. Two-pass to handle plates where heavy CR
     staining bleeds into the strict-red mask alongside the printed
     ring.
  2. `_find_ring_edges` snaps the precise inner / outer of every
     printed ring (red, gray, black) from the radial colour profile
     around the detected card centre. Spec ratios are starting
     points; observed colour transitions win.

`detect_geometry` orchestrates both, then runs `_detect_dish` to find
the petri-dish boundary independently from the printed-card centre.
"""
from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np

from .config import (
    DISH_BLACK_BAND_THRESH, HOUGH_DOWNSCALE, HOUGH_DP,
    N_SECTORS, RING_FRAC, SECTOR_0_THETA_DEG,
)


@dataclass
class Geometry:
    # Card / cutout: printed rings are concentric with this centre.
    card_cx: float
    card_cy: float
    r_card: float
    r_red_in: float
    r_red_out: float
    r_gray_in: float
    r_gray_out: float
    r_black_in: float
    r_black_out: float
    # Petri dish: where the agar and colonies live; can be offset
    # within the cutout's 5 mm clearance, so detected separately.
    dish_cx: float
    dish_cy: float
    r_dish: float


# ── Red-ring detection ────────────────────────────────────────

def _detect_red_ring(img: np.ndarray) -> tuple:
    """Return (card_cx, card_cy, r_red_inner, r_red_outer) by direct
    segmentation of the printed red ring. Robust to partial occlusion
    (Hough is voting-based) and to heavy CR staining bleeding into
    the strict-red mask (two-pass refinement keeps only the radially
    outermost component, which is geometrically always the printed
    ring since its outer diameter — 125 mm — exceeds the dish — 100
    mm).
    """
    h, w = img.shape[:2]
    small = cv2.resize(img, (int(w * HOUGH_DOWNSCALE), int(h * HOUGH_DOWNSCALE)),
                       interpolation=cv2.INTER_AREA)
    sh, sw = small.shape[:2]
    R = small[..., 0]
    G = small[..., 1]
    B = small[..., 2]
    eps = 1e-4
    # Saturation gate: printed red ink has R >> G and R >> B. Brown
    # CR-tinted agar has ratios ~1.5, so R/G ≥ 2 excludes the agar
    # while keeping the ring. A small absolute floor on R rejects
    # dark sensor noise that happens to satisfy the ratio gates.
    mask = ((R / (G + eps) >= 2.0)
            & (R / (B + eps) >= 2.0)
            & (R > 0.02))
    mask_u8 = (mask.astype(np.uint8)) * 255
    kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    mask_u8 = cv2.morphologyEx(mask_u8, cv2.MORPH_OPEN, kern)
    mask_u8 = cv2.morphologyEx(mask_u8, cv2.MORPH_CLOSE, kern)
    n_lab, lab, stats, _ = cv2.connectedComponentsWithStats(mask_u8 // 255, connectivity=8)
    if n_lab < 2:
        raise RuntimeError('red-ring detection: no saturated-red component found')
    # Pick the CC with the largest BOUNDING BOX (not raw area). On
    # heavily CR-stained plates the strict-red mask contains both
    # the printed ring and a large solid bleed inside the dish. The
    # bleed can have more pixels than the ring (the ring is hollow),
    # but its bbox is bounded by the dish diameter, while the ring's
    # bbox is the ring outer diameter — and the printed ring's outer
    # diameter always exceeds the dish (spec: 125 vs 100 mm). So
    # the largest-bbox component is always the printed ring.
    bbox_w = stats[1:, cv2.CC_STAT_WIDTH]
    bbox_h = stats[1:, cv2.CC_STAT_HEIGHT]
    bbox_extents = np.maximum(bbox_w, bbox_h)
    biggest = 1 + int(np.argmax(bbox_extents))
    ring_mask = (lab == biggest).astype(np.uint8) * 255
    if int(ring_mask.sum()) < 500 * 255:
        raise RuntimeError('red-ring detection: largest red component too small')

    # First Hough on the chosen CC: edges are clean and circular,
    # so this finds the ring centre robustly even when part of the
    # ring is occluded by the dish lid.
    short_s = min(sh, sw)
    circles = cv2.HoughCircles(
        ring_mask, cv2.HOUGH_GRADIENT,
        dp=HOUGH_DP, minDist=int(short_s * 0.5),
        param1=80, param2=20,
        minRadius=int(short_s * 0.18),
        maxRadius=int(short_s * 0.40),
    )
    if circles is None:
        raise RuntimeError('red-ring detection: Hough on red mask returned no circles')
    cx0, cy0, _ = circles[0][0]

    # Second-pass Hough on the radially-outermost portion of the
    # ring mask only. On heavily CR-stained plates the strict-red
    # mask contains both the printed ring AND a large bleed inside
    # the dish at smaller radii — a single Hough on the merged blob
    # fits a circle that compromises between them. The ring is
    # geometrically the OUTERMOST red feature, so keeping only the
    # top ~40 % of pixel radii from the rough centre lets the
    # second Hough lock onto the printed ring cleanly. Clean plates
    # lose nothing — their strict-red mask is mostly the ring.
    ys_all, xs_all = np.where(lab == biggest)
    rs0 = np.sqrt((xs_all - cx0) ** 2 + (ys_all - cy0) ** 2)
    r_cut = float(np.percentile(rs0, 60))
    keep = rs0 >= r_cut
    outer_mask = np.zeros_like(ring_mask)
    outer_mask[ys_all[keep], xs_all[keep]] = 255
    circles2 = cv2.HoughCircles(
        outer_mask, cv2.HOUGH_GRADIENT,
        dp=HOUGH_DP, minDist=int(short_s * 0.5),
        param1=80, param2=20,
        minRadius=int(short_s * 0.18),
        maxRadius=int(short_s * 0.40),
    )
    if circles2 is None:
        cx, cy, hough_r = cx0, cy0, circles[0][0][2]
    else:
        cx, cy, hough_r = circles2[0][0]

    # Estimate ring inner / outer from Hough's centerline radius
    # using the spec ratios (printed ring centerline = 115 mm; inner
    # = 105 mm, outer = 125 mm). Hough's voting-based centerline is
    # robust to red-mask leakage on heavily CR-stained plates;
    # `_find_ring_edges` snaps the precise edges from the radial
    # colour profile.
    SPEC_INNER_OVER_CENTER = 105.0 / 115.0
    SPEC_OUTER_OVER_CENTER = 125.0 / 115.0
    r_inner_s = float(hough_r) * SPEC_INNER_OVER_CENTER
    r_outer_s = float(hough_r) * SPEC_OUTER_OVER_CENTER
    return (float(cx) / HOUGH_DOWNSCALE,
            float(cy) / HOUGH_DOWNSCALE,
            r_inner_s / HOUGH_DOWNSCALE,
            r_outer_s / HOUGH_DOWNSCALE)


# ── Dish detection ────────────────────────────────────────────

def _detect_dish(img: np.ndarray, card_cx: float, card_cy: float,
                 r_card: float) -> tuple:
    """Return (dish_cx, dish_cy, r_dish) by finding the agar→black-band
    transition along radii from the card centre. Edge points are fit
    to a circle by algebraic least squares (Kasa method), with a
    single outlier-rejection refit pass.
    """
    h, w = img.shape[:2]
    small = cv2.resize(img, (int(w * HOUGH_DOWNSCALE), int(h * HOUGH_DOWNSCALE)),
                       interpolation=cv2.INTER_AREA)
    sh, sw = small.shape[:2]
    luminance = small.mean(axis=2)
    cx_s = card_cx * HOUGH_DOWNSCALE
    cy_s = card_cy * HOUGH_DOWNSCALE
    r_s = r_card * HOUGH_DOWNSCALE

    n_angles = 360
    r_lo = 0.60 * r_s
    r_hi = 1.00 * r_s
    rs = np.arange(r_lo, r_hi, 0.5)
    angles = np.linspace(0.0, 2 * np.pi, n_angles, endpoint=False)

    edge_xs, edge_ys = [], []
    for theta in angles:
        cos_t, sin_t = np.cos(theta), -np.sin(theta)
        xs = cx_s + rs * cos_t
        ys = cy_s + rs * sin_t
        valid = (xs >= 0) & (xs < sw) & (ys >= 0) & (ys < sh)
        if not valid.any():
            continue
        xi = xs[valid].astype(int)
        yi = ys[valid].astype(int)
        lums = luminance[yi, xi]
        lums_smooth = np.convolve(lums, np.ones(5) / 5, mode='same')
        # Scan from r_hi inward: agar is variably dark and full of
        # below-threshold patches, so an inside-out scan would
        # trigger inside the dish. Going outside-in we first cross
        # the uniformly very-dark black band, and the first sample
        # rising above the threshold is the agar (dish lid edge).
        n = len(lums_smooth)
        above = lums_smooth > DISH_BLACK_BAND_THRESH
        edge_idx = None
        for j in range(n - 1, -1, -1):
            if above[j]:
                # Require sustained agar (3+ consecutive samples) so
                # a single bright noise pixel doesn't trigger.
                if j >= 2 and above[j - 2:j + 1].all():
                    edge_idx = j
                    break
        if edge_idx is None:
            continue
        edge_xs.append(float(xi[edge_idx]))
        edge_ys.append(float(yi[edge_idx]))

    if len(edge_xs) < 30:
        raise RuntimeError(
            f'dish detection: only {len(edge_xs)} ray edges found '
            '(image too dim or no clear agar/black-band transition)'
        )
    xs_a = np.array(edge_xs)
    ys_a = np.array(edge_ys)

    A = np.column_stack([2 * xs_a, 2 * ys_a, np.ones_like(xs_a)])
    b = xs_a ** 2 + ys_a ** 2
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    dx, dy, k = sol
    r_fit = float(np.sqrt(k + dx * dx + dy * dy))
    resid = np.abs(np.sqrt((xs_a - dx) ** 2 + (ys_a - dy) ** 2) - r_fit)
    keep = resid <= 3.0
    if keep.sum() >= 30:
        xs_a = xs_a[keep]
        ys_a = ys_a[keep]
        A = np.column_stack([2 * xs_a, 2 * ys_a, np.ones_like(xs_a)])
        b = xs_a ** 2 + ys_a ** 2
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        dx, dy, k = sol
        r_fit = float(np.sqrt(k + dx * dx + dy * dy))

    return (float(dx) / HOUGH_DOWNSCALE,
            float(dy) / HOUGH_DOWNSCALE,
            r_fit / HOUGH_DOWNSCALE)


# ── Radial profile + ring-edge snap ───────────────────────────

def _radial_profile(img: np.ndarray, cx: float, cy: float,
                    r_max: int) -> np.ndarray:
    """Mean RGB per integer-radius bin from (cx, cy) out to r_max."""
    h, w = img.shape[:2]
    yy, xx = np.ogrid[:h, :w]
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    r_int = np.clip(r.astype(np.int32), 0, r_max).ravel()
    counts = np.bincount(r_int, minlength=r_max + 1).astype(np.float32)
    counts_safe = np.maximum(counts, 1.0)
    profile = np.empty((r_max + 1, 3), dtype=np.float32)
    for ch in range(3):
        sums = np.bincount(r_int,
                           weights=img[..., ch].astype(np.float32).ravel(),
                           minlength=r_max + 1)
        profile[:, ch] = sums / counts_safe
    profile[counts == 0] = 0.0
    return profile


# Predicate threshold for the red ring's snap pass — kept here as a
# named constant so the failure-diagnostic plotter can reference the
# same value without drift. Permissive (R-G > 0.020) to keep faint or
# heavily-CR-stained rings detectable; the longest-run selector below
# picks the printed ring even when other in-window radii also fire.
RED_PRED_THRESH = 0.020


def _find_ring_edges(profile: np.ndarray, r_card: float) -> dict:
    """Locate each printed ring's inner / outer edges from the radial
    profile around the card centre, using the spec ratios as starting
    points and snapping to the actual colour transitions.
    """
    n = profile.shape[0]
    R = profile[:, 0]
    G = profile[:, 1]
    B = profile[:, 2]
    lum = profile.mean(axis=1)

    edges = {}
    # Spec windows (fraction of r_card); widened by ±0.12 to absorb
    # printing / cropping variability around the nominal positions.
    spec = {
        'red':   (RING_FRAC['red'][0],   RING_FRAC['red'][1]),
        'gray':  (RING_FRAC['gray'][0],  RING_FRAC['gray'][1]),
        'black': (RING_FRAC['black'][0], RING_FRAC['black'][1]),
    }
    pred = {
        'red':   (R - G > RED_PRED_THRESH),
        'gray':  (np.abs(R - G) < 0.012) & (np.abs(G - B) < 0.012)
                 & (lum > 0.030) & (lum < 0.075),
        'black': (np.abs(R - G) < 0.012) & (np.abs(G - B) < 0.012)
                 & (lum < 0.025),
    }
    for name, (f_lo, f_hi) in spec.items():
        lo_search = int(max(0, (f_lo - 0.12) * r_card))
        hi_search = int(min(n - 1, (f_hi + 0.12) * r_card))
        if lo_search >= hi_search:
            raise RuntimeError(f'ring "{name}" search window empty')
        flag = pred[name][lo_search:hi_search + 1]
        if not flag.any():
            raise RuntimeError(
                f'ring "{name}" not found in radial profile '
                f'(searched {lo_search}-{hi_search} px)'
            )
        idx = np.where(flag)[0]
        breaks = np.where(np.diff(idx) > 1)[0]
        if breaks.size:
            starts = np.r_[0, breaks + 1]
            ends = np.r_[breaks, idx.size - 1]
            run_lens = ends - starts
            longest = int(np.argmax(run_lens))
            run_lo = idx[starts[longest]]
            run_hi = idx[ends[longest]]
        else:
            run_lo = idx[0]
            run_hi = idx[-1]
        edges[f'{name}_in'] = float(lo_search + run_lo)
        edges[f'{name}_out'] = float(lo_search + run_hi)
    return edges


def detect_geometry(img: np.ndarray) -> Geometry:
    """Detect card centre via the red ring, snap each printed ring's
    true inner / outer from the radial colour profile, and detect the
    petri dish independently. Ring positions come from observed
    colour transitions rather than spec × r_card, so small printing
    or lens variations are absorbed.
    """
    card_cx, card_cy, r_red_inner, r_red_outer = _detect_red_ring(img)
    r_card_seed = 0.5 * (r_red_inner / RING_FRAC['red'][0]
                         + r_red_outer / RING_FRAC['red'][1])
    h, w = img.shape[:2]
    r_max_profile = int(min(h, w) * 0.95)
    profile = _radial_profile(img, card_cx, card_cy, r_max_profile)
    edges = _find_ring_edges(profile, r_card_seed)
    r_card = 0.5 * (edges['red_in'] + edges['red_out']) / 1.15
    dish_cx, dish_cy, r_dish = _detect_dish(img, card_cx, card_cy, r_card)
    return Geometry(
        card_cx=card_cx, card_cy=card_cy, r_card=r_card,
        r_red_in=edges['red_in'], r_red_out=edges['red_out'],
        r_gray_in=edges['gray_in'], r_gray_out=edges['gray_out'],
        r_black_in=edges['black_in'], r_black_out=edges['black_out'],
        dish_cx=dish_cx, dish_cy=dish_cy, r_dish=r_dish,
    )


# ── Polar grids + sector mask ─────────────────────────────────

def _angle_grid(cx: float, cy: float, h: int, w: int) -> np.ndarray:
    """Per-pixel angle in [0, 360); 0° = right, 90° = top (image y flipped)."""
    yy, xx = np.ogrid[:h, :w]
    theta = np.degrees(np.arctan2(-(yy - cy), xx - cx))
    return np.where(theta < 0, theta + 360.0, theta).astype(np.float32)


def _radius_grid(cx: float, cy: float, h: int, w: int) -> np.ndarray:
    yy, xx = np.ogrid[:h, :w]
    return np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2).astype(np.float32)


def _sector_mask(theta_deg: np.ndarray, sector: int) -> np.ndarray:
    span = 360.0 / N_SECTORS
    half = span / 2.0
    center = (SECTOR_0_THETA_DEG + sector * span) % 360.0
    lo = (center - half) % 360.0
    hi = (center + half) % 360.0
    if lo < hi:
        return (theta_deg >= lo) & (theta_deg < hi)
    return (theta_deg >= lo) | (theta_deg < hi)
