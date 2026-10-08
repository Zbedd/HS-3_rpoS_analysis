"""Per-spot quantification: the OD-chromatic metric.

The annulus geometry (band offset OUTSIDE the colony mask via
`cv2.dilate`, 1.30–1.50·r_eff with floors) is the local reference for
`OD_chrom = OD_G − OD_R`, the assay's only metric. It is emitted on
every SpotMetrics record and is what the driver plots.
"""
from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np

from .config import (
    ANNULUS_INNER_FRAC, ANNULUS_MIN_GAP_PX, ANNULUS_MIN_THICK_PX,
    ANNULUS_OUTER_FRAC, EPS, PLATE_INTERIOR_FRAC,
)
from .geometry import Geometry, _radius_grid


# ── Spot metrics dataclass ────────────────────────────────────

@dataclass
class SpotMetrics:
    sample: str
    sector: int
    cx: float
    cy: float
    r_eff: float
    area_px: int
    bg_area_px: int
    # ── OD-chromatic (raw-domain log-ratios on segmented core +
    # annulus reference; differences-of-logs cancel exposure / WB).
    core_R: float
    core_G: float
    ref_R: float
    ref_G: float
    OD_R: float
    OD_G: float
    OD_chrom: float


# ── Mask + annulus helpers ────────────────────────────────────

def _mask_from_contour(contour, h: int, w: int) -> np.ndarray:
    mask = np.zeros((h, w), dtype=np.uint8)
    pts = np.asarray(contour, dtype=np.int32)
    if pts.size == 0:
        return mask.astype(bool)
    cv2.fillPoly(mask, [pts], 255)
    return mask.astype(bool)


def _annulus_offsets(r_eff: float) -> tuple:
    """(gap_px, thick_px) for the annulus around a colony of effective
    radius `r_eff`."""
    gap_px = int(max(np.ceil((ANNULUS_INNER_FRAC - 1.0) * r_eff),
                     ANNULUS_MIN_GAP_PX))
    thick_px = int(max(np.ceil((ANNULUS_OUTER_FRAC - ANNULUS_INNER_FRAC) * r_eff),
                       ANNULUS_MIN_THICK_PX))
    return gap_px, thick_px


def annulus_for_mask(mask: np.ndarray, r_eff: float) -> np.ndarray:
    """Annulus that follows the colony's mask shape via `cv2.dilate`.
    Returns a boolean mask of the band at [`ANNULUS_INNER_FRAC`,
    `ANNULUS_OUTER_FRAC`]·r_eff outside the colony."""
    gap_px, thick_px = _annulus_offsets(r_eff)
    outer_off = gap_px + thick_px
    m_u8 = mask.astype(np.uint8)
    kern_in = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (2 * gap_px + 1, 2 * gap_px + 1))
    kern_out = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE, (2 * outer_off + 1, 2 * outer_off + 1))
    ring = cv2.dilate(m_u8, kern_out).astype(bool) \
           & ~cv2.dilate(m_u8, kern_in).astype(bool)
    return ring & ~mask


# ── Per-spot quantification (OD-chromatic) ────────────────────

def quantify_spots(img_raw: np.ndarray, geom: Geometry,
                   spots: list) -> list:
    """Per-spot **OD-chromatic** metric on the raw linear-light image:

    `OD_chrom = log(ref_G/core_G) - log(ref_R/core_R)`, with `core_*` =
    mean over the segmented mask and `ref_*` = median over the same
    annulus. Differences-of-logs cancel exposure / white-balance
    scaling, so no per-angle calibration is needed.
    """
    h, w = img_raw.shape[:2]
    rgrid = _radius_grid(geom.dish_cx, geom.dish_cy, h, w)
    interior = rgrid < geom.r_dish * PLATE_INTERIOR_FRAC

    masks = [_mask_from_contour(s.contour, h, w) for s in spots]
    if len(masks) > 1:
        union_others = [
            np.logical_or.reduce([m for j, m in enumerate(masks) if j != i])
            for i in range(len(masks))
        ]
    else:
        union_others = [np.zeros((h, w), dtype=bool)] * len(masks)

    out = []
    for spot, mask, others in zip(spots, masks, union_others):
        ring = annulus_for_mask(mask, spot.r_eff)
        annulus = ring & interior & ~others
        if not mask.any():
            raise RuntimeError(f'core mask empty for sector {spot.sector}')
        if not annulus.any():
            raise RuntimeError(f'annulus empty for sector {spot.sector}')

        # ── OD-chromatic (raw-domain log-ratios). Computed on the raw
        # linear image, NOT a calibrated one: a black-anchored,
        # unclipped calibration sends agar background and the green
        # channel ≤ 0, which is undefined for a log reflectance ratio.
        # The raw image is always positive, and the adjacent ref/core
        # ratio already cancels local illumination and white balance.
        core_R = max(float(img_raw[..., 0][mask].mean()), EPS)
        core_G = max(float(img_raw[..., 1][mask].mean()), EPS)
        ref_R = max(float(np.median(img_raw[..., 0][annulus])), EPS)
        ref_G = max(float(np.median(img_raw[..., 1][annulus])), EPS)
        OD_R = float(np.log(ref_R / core_R))
        OD_G = float(np.log(ref_G / core_G))
        OD_chrom = OD_G - OD_R

        out.append(SpotMetrics(
            sample=spot.sample, sector=spot.sector,
            cx=spot.cx, cy=spot.cy,
            r_eff=spot.r_eff, area_px=spot.area_px,
            bg_area_px=int(annulus.sum()),
            core_R=core_R, core_G=core_G,
            ref_R=ref_R, ref_G=ref_G,
            OD_R=OD_R, OD_G=OD_G, OD_chrom=OD_chrom,
        ))
    return out
