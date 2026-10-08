"""I/O helpers: DNG decode, PNG save, JSON cache, CSV r/w.

Functions here either consume / produce bytes on disk or transform
between dict-of-strings (CSV row) and dict-of-typed-values (Python row)
representations. No metric calculation.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Optional

import cv2
import numpy as np
import rawpy

from .config import (
    CACHE_DIR, DATA_DIR, NAME_RE, OUT_DIR, QUANT_CSV, RAWPY_KWARGS,
    SUMMARY_CSV,
)


# ── Image I/O ──────────────────────────────────────────────────

def load_dng_linear(path: Path) -> np.ndarray:
    """DNG decoded to linear-light RGB float32 in [0, 1]."""
    with rawpy.imread(str(path)) as raw:
        rgb16 = raw.postprocess(**RAWPY_KWARGS)
    return rgb16.astype(np.float32) / 65535.0


def save_png16(img_float: np.ndarray, out_path: Path) -> None:
    arr = (np.clip(img_float, 0, 1) * 65535.0 + 0.5).astype(np.uint16)
    bgr = cv2.cvtColor(arr, cv2.COLOR_RGB2BGR)
    if not cv2.imwrite(str(out_path), bgr):
        raise IOError(f'cv2.imwrite failed for {out_path}')


def parse_meta(stem: str) -> dict:
    m = NAME_RE.match(stem)
    if m is None:
        raise ValueError(f'unrecognized filename stem: {stem!r}')
    return {
        'time': m.group('time'),
        'media': m.group('media'),
        'plate_index': int(m.group('idx')),
    }


# ── Detection cache ────────────────────────────────────────────

def cache_path(stem: str) -> Path:
    return CACHE_DIR / f'{stem}.json'


def load_cache(stem: str) -> Optional[dict]:
    p = cache_path(stem)
    return json.loads(p.read_text()) if p.exists() else None


def save_cache(stem: str, payload: dict) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path(stem).write_text(json.dumps(payload, indent=2))


def cached_r_eff(stem: str, sample: str) -> Optional[float]:
    """`r_eff` of one cached colony, or None when it isn't cached."""
    cache = load_cache(stem)
    if cache is None:
        return None
    spot = next((s for s in cache['spots'] if s['sample'] == sample), None)
    return None if spot is None else float(spot['r_eff'])


def load_colony_crops(stem: str, samples, half_px: int) -> dict:
    """Square crops of the named colonies from the plate's raw linear
    image, `half_px` either side of each cached centre. Returns
    `sample -> linear RGB float array in [0, 1]`; a sample is absent
    when the DNG or its cached spot is missing. `half_px` is a plain
    pixel count, so every crop shares one scale.

    Raw, not the normalized PNG: the black-anchored calibration drives
    the green channel below zero and `save_png16` clips it there, so a
    normalized crop carries none of the green that `OD_chrom` is half
    built from. `quantify_spots` measures on the raw image for the same
    reason.
    """
    dng = DATA_DIR / f'{stem}.DNG'
    cache = load_cache(stem)
    if not dng.exists() or cache is None:
        return {}
    plate = load_dng_linear(dng)
    h, w = plate.shape[:2]
    by_sample = {s['sample']: s for s in cache['spots']}
    crops = {}
    for name in samples:
        spot = by_sample.get(name)
        if spot is None:
            continue
        cx, cy = int(round(spot['cx'])), int(round(spot['cy']))
        x0, x1 = max(0, cx - half_px), min(w, cx + half_px)
        y0, y1 = max(0, cy - half_px), min(h, cy + half_px)
        if x1 - x0 < 2 or y1 - y0 < 2:
            continue
        crops[name] = plate[y0:y1, x0:x1].copy()
    return crops


# ── Quantification CSV ─────────────────────────────────────────

QUANT_FIELDS = [
    'image', 'time', 'media', 'plate_index', 'sample', 'sector',
    'cx', 'cy', 'r_eff', 'area_px', 'bg_area_px',
    # OD-chromatic
    'core_R', 'core_G', 'ref_R', 'ref_G',
    'OD_R', 'OD_G', 'OD_chrom',
]

SUMMARY_FIELDS = [
    'image', 'time', 'media', 'plate_index', 'sample',
    'OD_chrom',
]

_QUANT_INT_KEYS = ('plate_index', 'sector', 'area_px', 'bg_area_px')
_QUANT_FLOAT_KEYS = (
    'cx', 'cy', 'r_eff',
    'core_R', 'core_G', 'ref_R', 'ref_G',
    'OD_R', 'OD_G', 'OD_chrom',
)


def write_csvs(rows: list) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with QUANT_CSV.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=QUANT_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in QUANT_FIELDS})
    with SUMMARY_CSV.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in SUMMARY_FIELDS})
    print(f'wrote {QUANT_CSV} ({len(rows)} rows)')
    print(f'wrote {SUMMARY_CSV}')


def read_summary() -> list:
    if not SUMMARY_CSV.exists():
        sys.exit(f'summary CSV not found at {SUMMARY_CSV}; '
                 'run without --plots-only first')
    float_cols = ('OD_chrom',)
    with SUMMARY_CSV.open(newline='') as f:
        rows = []
        for r in csv.DictReader(f):
            r['plate_index'] = int(r['plate_index'])
            for k in float_cols:
                if k in r and r[k] != '':
                    r[k] = float(r[k])
            rows.append(r)
        return rows


def _cast_quant_row(r: dict) -> dict:
    """Cast a CSV-loaded quantification row back to native types so it
    can be passed straight to plot_*."""
    out = dict(r)
    for k in _QUANT_INT_KEYS:
        out[k] = int(out[k])
    for k in _QUANT_FLOAT_KEYS:
        if k in out and out[k] != '':
            out[k] = float(out[k])
    return out


def load_existing_quant() -> dict:
    """Map of `image stem -> list of typed rows` from a prior run's
    quantification.csv. Returns {} if the CSV doesn't exist yet."""
    if not QUANT_CSV.exists():
        return {}
    by_image = {}
    with QUANT_CSV.open(newline='') as f:
        for r in csv.DictReader(f):
            by_image.setdefault(r['image'], []).append(_cast_quant_row(r))
    return by_image
