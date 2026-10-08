"""Congo Red pipeline configuration.

All module-level constants used across the pipeline live here, so the
science modules don't carry their own copies. `quantify_spots` computes
`OD_chrom`, the assay's only metric.
"""
from __future__ import annotations

import re

import numpy as np

import paths
import viz as fs


# ── Paths ──────────────────────────────────────────────────────

DATA_DIR = paths.data('CR')
OUT_DIR = paths.outputs('CR')
NORM_DIR = OUT_DIR / 'normalized'
ANNOT_DIR = OUT_DIR / 'annotated'
CACHE_DIR = OUT_DIR / 'cache'
FAILURES_DIR = OUT_DIR / 'failures'
PLOTS_DIR = OUT_DIR / 'plots'

QUANT_CSV = OUT_DIR / 'quantification.csv'
SUMMARY_CSV = OUT_DIR / 'cr_summary.csv'

# ── Plot toggles ──────────────────────────────────────────────

# Colony crops above each panel of the effect-estimation figure. Needs the
# raw DNGs and the detection cache; the figure drops the strip when
# either is absent.
SHOW_COLONY_STRIP = True


# ── Decoding ───────────────────────────────────────────────────

# Linear-light DNG decode.
RAWPY_KWARGS = dict(
    gamma=(1, 1),
    no_auto_bright=True,
    output_bps=16,
    use_camera_wb=True,
)


# ── Reference-card geometry (mm) ───────────────────────────────

DISH_MM = 100.0
RING_DIAM_MM = {  # (inner_diameter, outer_diameter)
    'red':   (105.0, 125.0),
    'gray':  (130.0, 150.0),
    'black': (155.0, 175.0),
}
RING_FRAC = {name: (i / DISH_MM, o / DISH_MM)
             for name, (i, o) in RING_DIAM_MM.items()}


def _srgb_to_linear_scalar(srgb_uint8: int) -> float:
    v = srgb_uint8 / 255.0
    return v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4


# Calibration anchor truths (linear light, [0, 1]).
RING_TRUTH = {
    'black': np.zeros(3, dtype=np.float32),
    'gray':  np.full(3, _srgb_to_linear_scalar(128), dtype=np.float32),
}


# ── Sectoring ──────────────────────────────────────────────────

# Six 60° wedges; sector 0 is the blank (top of the plate, θ = 90°);
# sector indices increase CCW. Misorientated plates must be reshot or
# rotated externally to match this layout.
N_SECTORS = 6
SECTOR_0_THETA_DEG = 90.0
# ⚠ DATA EXCLUSION — sector 1 (`WT`, the unmarked parental isolate, drawn
# as Parental) is contaminated, do NOT use in analysis. Its R2A plates were
# bacterially contaminated, so the condition is unreliable and must be
# excluded from every figure / statistic / comparison. It is still detected
# + quantified and kept in the CSVs for provenance only. (See the header
# note in cr_pipeline.py.)
SECTOR_SAMPLE = {
    0: 'blank',
    1: 'WT',
    2: 'kanR',
    3: 'rpoS',
    4: 'PA14',
    5: 'PA14_rpoS',
}
# Every non-blank sector. This is the *quantification* list: each of
# these is detected, measured and written to the CSVs.
SAMPLES = [SECTOR_SAMPLE[i] for i in range(1, N_SECTORS)]
BLANK_SECTOR = 0

# Samples kept in the CSVs for provenance but barred from every figure
# and statistic. The parental isolate had contamination. Anything listed
# here must never reach a plot.
EXCLUDED_SAMPLES = ('WT',)

# The *analysis* list — what figures and statistics may use. Derived
# from SAMPLES so a new sector is picked up automatically, and from
# EXCLUDED_SAMPLES so an exclusion only has to be declared once.
PLOT_SAMPLES = [s for s in SAMPLES if s not in EXCLUDED_SAMPLES]


# ── Geometry detection ────────────────────────────────────────

HOUGH_DP = 1.2
HOUGH_DOWNSCALE = 0.25  # downsample factor for the Hough passes

# Dish edge: the cutout-vs-dish-lid gap shows the black table behind
# the card (uniformly dark). Each radial ray crosses this band only
# when it leaves the agar; the first sustained luminance dip marks the
# dish lid edge. Set comfortably below typical agar luminance (~0.02+)
# and above sensor noise (~0.005).
DISH_BLACK_BAND_THRESH = 0.015

# Calibration QC: reject only inverted gaps (gray darker than black,
# meaning detection sampled the wrong rings); a small positive gap is
# fine — it just gives a higher slope and amplifies sensor noise
# proportionally. Empirical gaps run 0.008–0.030.
CALIB_MIN_GAP = 0.005


# ── Spot detection ─────────────────────────────────────────────

SPOT_R_INNER_FRAC = 0.18      # search annulus inner edge (frac of r_dish)
SPOT_R_OUTER_FRAC = 0.72
SPOT_DIST_PERCENTILE = 96.0   # default threshold percentile of in-wedge colour distance
SPOT_MORPH_KSIZE = 5          # close→open kernel size (px) for the initial cleanup
SPOT_SMOOTH_KSIZE = 7         # post-CC CLOSE kernel — smooths concavities; CLOSE only
                              # (no OPEN) so thin parts of small colonies aren't eroded
SPOT_MIN_AREA_FRAC = 0.0005   # of dish area (~1700 px) — below = noise, skip
SPOT_MAX_AREA_FRAC = 0.05
SPOT_CENTER_ANGLE_TOL_DEG = 25.0  # pick candidate within this of wedge centerline

# Quality floor: ~70% of the typical post-detection median area. Spots
# below this are kept but flagged `below_floor=True` in the cache so
# downstream code can warn or filter. Real-but-small colonies are
# tagged rather than crashing the run.
SPOT_MIN_AREA_PX = 7500

# Auto-override threshold for the crescent failure mode: a colony that fills
# enough of its wedge biases the per-wedge median toward its own colour, so
# the default segments only the rim. Solidity (contour area / convex-hull
# area) separates these — solid blobs ≥ 0.95, elongated colonies ≥ 0.90,
# crescents ≤ 0.80. Spots below the floor are auto-retried; highest solidity
# wins.
SPOT_SOLIDITY_FLOOR = 0.85
# Spots with area < SPOT_AREA_RATIO_FLOOR × plate-median are also
# auto-retried — covers cases where rim-only detection produces a
# convex but much-smaller-than-peers shape.
SPOT_AREA_RATIO_FLOOR = 0.55

# Sanity cap on Otsu (used by `blank_otsu`): colonies cover ~5–10 % of
# a wedge in this layout; if Otsu would select more than 15 % the
# histogram lacks a clean bimodal split, so fall back to a percentile.
SPOT_OTSU_MAX_FRAC = 0.15
SPOT_OTSU_FALLBACK_PERCENTILE = 96.0


# ── Annulus geometry (post-detection sampling) ─────────────────

# Implemented as a band offset OUTSIDE the colony's actual mask via
# `cv2.dilate`, so it follows non-circular shapes. Floors keep the band
# usable on small / irregular colonies. Used by the OD-chromatic
# quantification and the annotation overlays.
ANNULUS_INNER_FRAC = 1.30
ANNULUS_OUTER_FRAC = 1.50
ANNULUS_MIN_GAP_PX = 10
ANNULUS_MIN_THICK_PX = 14

# Dish-interior fraction. The annulus is intersected with this disc to
# avoid sampling the dish-edge / lid region.
PLATE_INTERIOR_FRAC = 0.92


# ── Per-spot manual overrides ──────────────────────────────────

# Each entry runs an alternative detection on just the listed
# (stem, sample) pair after the default pass; everything else on the
# plate keeps the default settings.
#
# Available `method` values:
#   "blank_otsu"  — blank-sector (clean agar) median + Otsu threshold with a
#                   sanity cap. For the crescent mode described above.
#   "raw_p96"     — default logic (per-wedge median + p96) on the raw
#                   linear-light image. For plates where calibration squashes
#                   the dish to near-zero negatives and detection finds noise.
#   "convex_hull" — `raw_p96` plus the convex hull of the contour. For
#                   colonies whose rim is the only high-contrast feature, so
#                   every pixel-level method reads a crescent.
#
# New entries need `--redetect` to take effect. `detect_colonies` also
# auto-retries alternatives for spots with low solidity / area; entries here
# force a method when that auto-pick is wrong.
TARGETED_ADJUSTMENTS = {
    '48h_LB_2': {
        'PA14_rpoS': {'method': 'blank_otsu'},
    },
    '72h_NaCl-_5': {
        'kanR': {'method': 'convex_hull'},
    },
    '72h_NaCl-_12': {
        'rpoS': {'method': 'convex_hull'},
    },
}

# Per-(stem, sample) entries that should be skipped entirely — colonies
# that physically didn't grow / weren't spotted, where any "detection"
# would be noise. Treated as NA: detection is skipped, the spot is
# omitted from all CSVs and downstream metric computations, and the
# annotation overlay marks the sector with an "NA" label so the
# missing entry is visible (distinct from "detection silently
# failed"). Any sector not listed here goes through normal detection.
#
#   stem -> set of sample names to skip
NA_SPOTS = {
    # Same physical plate at two timepoints — the PA14_rpoS sector
    # never produced a visible colony at either reading.
    '48h_R2A_4': {'PA14_rpoS'},
    '72h_R2A_9': {'PA14_rpoS'},
}


# ── Plotting / aggregation ─────────────────────────────────────

TIMES = ['48h', '72h']
MEDIA_ORDER = ['LB', 'NaCl-', 'R2A']
MEDIA_DISPLAY = {'LB': 'Miller LB', 'NaCl-': 'NaCl-free LB', 'R2A': 'R2A'}
# Media quantified into the CSVs but barred from every figure — NaCl-free LB
# was run for this assay alone, so it has no counterpart in the other assays.
EXCLUDED_MEDIA = ('NaCl-',)
PLOT_MEDIA = [m for m in MEDIA_ORDER if m not in EXCLUDED_MEDIA]
# Strain naming is the whole repository's, so `viz` owns it — see STRAIN
# NAMING in `viz/style.py`.
RPOS_MINUS = fs.RPOS_MINUS
RNG_SEED = 0
EPS = 1e-4

# The rpoS-disruption effect, one series per strain background. The parental
# isolate is absent (contaminated — see the sector map above); each series
# contrasts the rpoS mutant against its own marker/parent control. No
# `color`: a difference takes a contrast hue by position, held apart from the
# strain colours it is a difference of.
RPOS_EFFECT = [
    {'label': 'HS-3', 'numerator': 'rpoS', 'denominator': 'kanR'},
    {'label': 'PA14', 'numerator': 'PA14_rpoS', 'denominator': 'PA14'},
]

NAME_RE = re.compile(r'^(?P<time>48h|72h)_(?P<media>LB|NaCl-|R2A)_(?P<idx>\d+)$')
