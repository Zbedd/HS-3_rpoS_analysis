"""Annotation overlay + per-plate failure diagnostic.

`annotate` draws the dish, printed rings, sector boundaries, blank
marker, per-spot colony outline, annulus bands, sample / metric
labels, and any NA-sector markers onto the calibrated image.

`save_failure_diagnostic` builds a multi-panel PNG showing whatever
geometry was successfully extracted, the strict-red mask + chosen
component, and (if available) the radial profile around the detected
card centre. Designed so the failure mode is visible at a glance
when re-running diagnoses in `outputs/CR/failures/`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import cv2
import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

from .config import (
    BLANK_SECTOR, FAILURES_DIR, HOUGH_DOWNSCALE, N_SECTORS,
    RING_FRAC, SECTOR_0_THETA_DEG, SECTOR_SAMPLE,
)
from .geometry import Geometry, _detect_red_ring, _radial_profile
from .geometry import RED_PRED_THRESH
from .io import load_dng_linear
from .quantify import _annulus_offsets, _mask_from_contour


# ── Display helpers ───────────────────────────────────────────

def _img_for_display(calibrated: np.ndarray) -> np.ndarray:
    """Linear → gamma 2.2 + 99th-percentile stretch (for viewing only)."""
    g = np.clip(calibrated, 0, 1) ** (1.0 / 2.2)
    p99 = float(np.percentile(g, 99))
    if p99 > 1e-6:
        g = np.clip(g / p99, 0, 1)
    return (g * 255.0 + 0.5).astype(np.uint8)


def _label(bgr, text, pos, color, scale, stroke, fill):
    cv2.putText(bgr, text, pos, cv2.FONT_HERSHEY_SIMPLEX, scale, (0, 0, 0),
                stroke, cv2.LINE_AA)
    cv2.putText(bgr, text, pos, cv2.FONT_HERSHEY_SIMPLEX, scale, color,
                fill, cv2.LINE_AA)


# ── Annotation overlay ───────────────────────────────────────

def annotate(calibrated: np.ndarray, geom: Geometry,
             spots: list, metrics: list, stem: str,
             na_samples: Optional[set] = None) -> np.ndarray:
    h, w = calibrated.shape[:2]
    disp = _img_for_display(calibrated)
    bgr = cv2.cvtColor(disp, cv2.COLOR_RGB2BGR)
    card_cx_i, card_cy_i = int(geom.card_cx), int(geom.card_cy)
    dish_cx_i, dish_cy_i = int(geom.dish_cx), int(geom.dish_cy)

    # Dish (centered on dish_cx/cy).
    cv2.circle(bgr, (dish_cx_i, dish_cy_i), int(geom.r_dish),
               (255, 255, 0), 3, cv2.LINE_AA)
    cv2.drawMarker(bgr, (dish_cx_i, dish_cy_i), (255, 255, 0),
                   cv2.MARKER_CROSS, 30, 2)

    # Rings (concentric with the printed card centre).
    ring_bgr = {'red': (40, 40, 220), 'gray': (160, 160, 160), 'black': (60, 60, 60)}
    for name, color in ring_bgr.items():
        ri = int(getattr(geom, f'r_{name}_in'))
        ro = int(getattr(geom, f'r_{name}_out'))
        cv2.circle(bgr, (card_cx_i, card_cy_i), ri, color, 2, cv2.LINE_AA)
        cv2.circle(bgr, (card_cx_i, card_cy_i), ro, color, 2, cv2.LINE_AA)

    # Sector boundaries (rays from dish centre — sectors are around
    # the dish, where samples were spotted).
    span = 360.0 / N_SECTORS
    for i in range(N_SECTORS):
        boundary = (SECTOR_0_THETA_DEG - span / 2.0 + i * span) % 360.0
        rad = np.radians(boundary)
        x2 = int(geom.dish_cx + geom.r_dish * np.cos(rad))
        y2 = int(geom.dish_cy - geom.r_dish * np.sin(rad))
        cv2.line(bgr, (dish_cx_i, dish_cy_i), (x2, y2),
                 (140, 140, 140), 1, cv2.LINE_AA)

    # Blank-sector marker.
    blank_theta = (SECTOR_0_THETA_DEG + BLANK_SECTOR * span) % 360.0
    blank_rad = np.radians(blank_theta)
    bx = int(geom.dish_cx + 0.55 * geom.r_dish * np.cos(blank_rad))
    by = int(geom.dish_cy - 0.55 * geom.r_dish * np.sin(blank_rad))
    _label(bgr, 'blank', (bx - 60, by), (200, 200, 200), 1.4, 8, 3)

    # Per-spot: colony outline, annulus rings, sample name + OD_chrom.
    for spot, m in zip(spots, metrics):
        pts = np.asarray(spot.contour, dtype=np.int32).reshape(-1, 1, 2)
        cv2.drawContours(bgr, [pts], -1, (0, 255, 0), 3, cv2.LINE_AA)

        # Annulus inner / outer dilation boundaries — same offsets
        # used by quantification.
        gap_px, thick_px = _annulus_offsets(spot.r_eff)
        outer_off = gap_px + thick_px
        m_u8 = _mask_from_contour(spot.contour, h, w).astype(np.uint8)
        for off in (gap_px, outer_off):
            kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,
                                             (2 * off + 1, 2 * off + 1))
            dil = cv2.dilate(m_u8, kern)
            ann_contours, _ = cv2.findContours(dil, cv2.RETR_EXTERNAL,
                                                cv2.CHAIN_APPROX_SIMPLE)
            cv2.drawContours(bgr, ann_contours, -1, (0, 255, 0), 1, cv2.LINE_AA)

        nx = int(spot.cx - spot.r_eff)
        ny = int(spot.cy - spot.r_eff - 50)
        _label(bgr, spot.sample, (nx, ny), (0, 255, 0), 1.8, 8, 4)
        _label(bgr, f'OD_chrom={m.OD_chrom:+.3f}',
               (nx, ny + 38), (0, 255, 0), 1.0, 5, 2)

    # NA-sector markers — draw a clear "NA" + sample name in the wedge
    # centre so an intentionally-skipped sector is visibly different
    # from a silent failure.
    if na_samples:
        for sample in na_samples:
            try:
                sec = next(s for s, n in SECTOR_SAMPLE.items() if n == sample)
            except StopIteration:
                continue
            sec_theta = (SECTOR_0_THETA_DEG + sec * span) % 360.0
            rad = np.radians(sec_theta)
            x = int(geom.dish_cx + 0.55 * geom.r_dish * np.cos(rad))
            y = int(geom.dish_cy - 0.55 * geom.r_dish * np.sin(rad))
            _label(bgr, sample, (x - 90, y - 20),
                   (160, 160, 255), 1.4, 8, 3)
            _label(bgr, 'NA', (x - 30, y + 35),
                   (160, 160, 255), 1.6, 8, 4)

    _label(bgr, stem, (15, 40), (255, 255, 255), 0.9, 4, 2)
    return bgr


# ── Failure diagnostic ───────────────────────────────────────

def save_failure_diagnostic(dng_path: Path, error_msg: str) -> None:
    """Render a multi-panel debug PNG for plates that threw during
    detection. Shows whatever geometry was successfully extracted plus
    the strict-red mask and (if available) the radial colour profile.
    Mirrors `_detect_red_ring`'s bbox-extent ranking so the chosen
    component is the one detection actually picked."""
    stem = dng_path.stem
    try:
        img = load_dng_linear(dng_path)
    except Exception as e:
        print(f'[FAIL-DIAG] {stem}: cannot load DNG ({e})')
        return

    h, w = img.shape[:2]

    card_cx = card_cy = r_red_in = r_red_out = None
    try:
        card_cx, card_cy, r_red_in, r_red_out = _detect_red_ring(img)
    except Exception:
        pass

    # Rebuild the strict-red mask + chosen CC for visualization.
    small = cv2.resize(img, (int(w * HOUGH_DOWNSCALE), int(h * HOUGH_DOWNSCALE)),
                       interpolation=cv2.INTER_AREA)
    R_s, G_s, B_s = small[..., 0], small[..., 1], small[..., 2]
    eps = 1e-4
    mask_raw = ((R_s / (G_s + eps) >= 2.0)
                & (R_s / (B_s + eps) >= 2.0)
                & (R_s > 0.02))
    mask_u8 = (mask_raw.astype(np.uint8)) * 255
    kern = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    mask_u8 = cv2.morphologyEx(mask_u8, cv2.MORPH_OPEN, kern)
    mask_u8 = cv2.morphologyEx(mask_u8, cv2.MORPH_CLOSE, kern)
    n_lab, lab, stats, _ = cv2.connectedComponentsWithStats(mask_u8 // 255, connectivity=8)
    largest_mask = np.zeros_like(mask_raw)
    if n_lab >= 2:
        bbox_w = stats[1:, cv2.CC_STAT_WIDTH]
        bbox_h = stats[1:, cv2.CC_STAT_HEIGHT]
        bbox_extents = np.maximum(bbox_w, bbox_h)
        biggest = 1 + int(np.argmax(bbox_extents))
        largest_mask = (lab == biggest)
    mask_viz = np.zeros((mask_raw.shape[0], mask_raw.shape[1], 3), dtype=np.float32)
    mask_viz[mask_raw] = [0.40, 0.40, 0.40]
    mask_viz[largest_mask] = [1.00, 0.20, 0.20]

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])

    ax_img = fig.add_subplot(gs[0, 0])
    disp = np.clip(img, 0, 1) ** (1.0 / 2.2)
    p99 = float(np.percentile(disp, 99))
    if p99 > 1e-6:
        disp = np.clip(disp / p99, 0, 1)
    ax_img.imshow(disp)
    if card_cx is not None:
        ax_img.plot(card_cx, card_cy, marker='x', color='red', ms=18, mew=3)
        for r in (r_red_in, r_red_out):
            ax_img.add_patch(plt.Circle((card_cx, card_cy), r,
                                         fill=False, color='red', lw=1.5))
        ax_img.set_title(f'red-ring detection: OK  '
                         f'cx={card_cx:.0f} cy={card_cy:.0f}  '
                         f'r_in={r_red_in:.0f} r_out={r_red_out:.0f}',
                         fontsize=10)
    else:
        ax_img.set_title('red-ring detection: FAILED', fontsize=10)
    ax_img.axis('off')

    ax_mask = fig.add_subplot(gs[0, 1])
    ax_mask.imshow(mask_viz)
    ax_mask.set_title(
        f'strict-red mask  raw={int(mask_raw.sum())}px  '
        f'largest_bbox_CC={int(largest_mask.sum())}px  '
        f'(gray = predicate hits, red = chosen component by bbox extent)',
        fontsize=10)
    ax_mask.axis('off')

    ax_prof = fig.add_subplot(gs[1, :])
    if card_cx is not None:
        r_max_profile = int(min(h, w) * 0.95)
        profile = _radial_profile(img, card_cx, card_cy, r_max_profile)
        rs = np.arange(profile.shape[0])
        ax_prof.plot(rs, profile[:, 0], color='red', lw=0.8, label='R')
        ax_prof.plot(rs, profile[:, 1], color='green', lw=0.8, label='G')
        ax_prof.plot(rs, profile[:, 2], color='blue', lw=0.8, label='B')
        ax_prof.plot(rs, profile[:, 0] - profile[:, 1],
                     color='purple', lw=1.4, label='R − G  (red predicate)')
        ax_prof.axhline(RED_PRED_THRESH, color='purple', ls='--', lw=0.8, alpha=0.6,
                        label=f'R−G > {RED_PRED_THRESH:.3f} (red threshold)')
        r_card_seed = 0.5 * (r_red_in / RING_FRAC['red'][0]
                             + r_red_out / RING_FRAC['red'][1])
        spec_color = {'red': 'red', 'gray': 'gray', 'black': 'black'}
        for name, (f_lo, f_hi) in RING_FRAC.items():
            ax_prof.axvspan(f_lo * r_card_seed, f_hi * r_card_seed,
                            alpha=0.10, color=spec_color[name],
                            label=f'spec {name} band')
        ax_prof.set_xlabel('radius from card centre (px)')
        ax_prof.set_ylabel('linear-light reflectance')
        ax_prof.set_title('radial profile from card centre')
        ax_prof.legend(loc='upper right', ncol=4, fontsize=8)
        ax_prof.set_xlim(0, min(r_max_profile, int(1.3 * r_card_seed)))
        ax_prof.grid(True, alpha=0.3)
    else:
        ax_prof.text(0.5, 0.5,
                     'radial profile unavailable: red-ring detection failed',
                     ha='center', va='center', transform=ax_prof.transAxes,
                     fontsize=14, fontstyle='italic')
        ax_prof.axis('off')

    fig.suptitle(f'{stem}  —  FAILURE\n{error_msg}',
                 fontsize=11, fontweight='bold', y=0.995)
    fig.tight_layout()

    FAILURES_DIR.mkdir(parents=True, exist_ok=True)
    out_path = FAILURES_DIR / f'{stem}.png'
    try:
        fig.savefig(out_path, dpi=120, bbox_inches='tight')
        print(f'[FAIL-DIAG] wrote {out_path}')
    except Exception as e:
        print(f'[FAIL-DIAG] could not save {out_path}: {e}')
    finally:
        plt.close(fig)
