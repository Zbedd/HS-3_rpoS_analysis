"""Congo Red assay pipeline — single driver script.

⚠ DATA EXCLUSION — WT is contaminated, do NOT use in analysis.
    The WT R2A plates were bacterially contaminated, which makes the
    WT (sector 1) condition unreliable across the board. WT is still
    decoded / detected / quantified and left in the CSVs for
    provenance, but it MUST be excluded from any figure, statistic, or
    comparison. Figures draw from `PLOT_SAMPLES`, which omits every name
    in `EXCLUDED_SAMPLES`; the lollipops contrast each rpoS knockout
    against its own control (HS-3: rpoS- − kanR; PA14: rpoS- − PA14).
    See the WT note in `assays.cr.config` beside the sector map.

Per-plate orchestration: decode → geometry → calibration → detection
→ quantification (OD-chromatic) → annotation.
Iterates the full batch under `data/CR/` and emits CSVs +
per-plate annotated PNGs + the manuscript figure.

`OD_chrom` is the assay's only metric. The science lives in the `cr`
package (`assays.cr`); this file holds CLI + I/O glue only.

CLI
---
    python scripts/cr_pipeline.py                 # full batch (uses cache)
    python scripts/cr_pipeline.py --redetect      # force re-detection
    python scripts/cr_pipeline.py --only STEM     # one image (no extension)
    python scripts/cr_pipeline.py --plots-only    # regenerate plots from CSVs
    python scripts/cr_pipeline.py --panels        # also write standalone panels
    python scripts/cr_pipeline.py --skip-existing # only plates without an
                                                  # annotated PNG (reshoots)
    python scripts/cr_pipeline.py --rerender      # re-render normalized +
                                                  # annotated PNGs regardless
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

import cv2
import numpy as np

from assays.cr.annotate import annotate, save_failure_diagnostic
from assays.cr.calibration import (
    apply_calibration, build_calibration, sample_ring_per_angle,
)
from assays.cr.config import (
    ANNOT_DIR, CACHE_DIR, DATA_DIR, FAILURES_DIR,
    NA_SPOTS, NORM_DIR, OUT_DIR, PLOTS_DIR, RPOS_EFFECT,
)
from assays.cr.detection import Spot, detect_colonies
from assays.cr.geometry import Geometry, detect_geometry
from assays.cr.io import (
    load_dng_linear, load_cache, load_existing_quant, parse_meta,
    read_summary, save_cache, save_png16, write_csvs,
)
from assays.cr import figure as cr_figure
from assays.cr.quantify import quantify_spots
from figures import cr_rpoS_effect


# ── Figure set ─────────────────────────────────────────────────

# Holds every per-dish difference at both timepoints (-0.247 .. +0.060), with
# headroom above for the significance marks.
CR_EFFECT_YLIM = (-0.27, 0.10)


# The manuscript uses 72 h. The same plates were imaged at 48 h, but image
# numbering does not establish plate identity across the two readings.
FIGURE_TIMES = ['72h']


def plot_all(rows: list, out_dir: Path, show_strip: bool = None,
             *, panels: bool = False) -> None:
    """The manuscript figure and, optionally, standalone diagnostic panels."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if panels:
        cr_figure.plot_colonies(rows, RPOS_EFFECT, out_dir / 'cr_colonies.png',
                                times=FIGURE_TIMES)
        cr_figure.plot_levels(rows, RPOS_EFFECT, out_dir / 'cr_levels.png',
                              times=FIGURE_TIMES)
        cr_figure.plot_effect(rows, RPOS_EFFECT, out_dir / 'cr_effect.png',
                              ylim=CR_EFFECT_YLIM, times=FIGURE_TIMES)
    cr_rpoS_effect.plot(rows, RPOS_EFFECT, out_dir / 'cr_rpoS_effect.png',
                        times=FIGURE_TIMES, show_strip=show_strip)

# ── Per-image driver ───────────────────────────────────────────

def process_image(dng_path: Path, redetect: bool,
                  rerender: bool = False) -> tuple:
    stem = dng_path.stem
    meta = parse_meta(stem)
    img = load_dng_linear(dng_path)
    cache = None if redetect else load_cache(stem)

    na_samples = set(NA_SPOTS.get(stem, set()))

    if cache is None:
        print(f'[detect] {stem}')
        geom = detect_geometry(img)
        black = sample_ring_per_angle(img, geom, 'black')
        gray = sample_ring_per_angle(img, geom, 'gray')
        calibration = build_calibration(black, gray)
        calibrated = apply_calibration(img, geom, calibration)
        spots = detect_colonies(calibrated, geom, stem=stem, img_raw=img)
        save_cache(stem, {
            'geometry': asdict(geom),
            'calibration': calibration.tolist(),
            'spots': [asdict(s) for s in spots],
            'na_samples': sorted(na_samples),
        })
    else:
        print(f'[cache]  {stem}')
        geom = Geometry(**cache['geometry'])
        calibration = np.asarray(cache['calibration'], dtype=np.float32)
        spots = [Spot(**s) for s in cache['spots']]
        # Caches written before the NA-spot config have no na_samples; fall
        # back to the cached value only when the config supplies none.
        if not na_samples:
            na_samples = set(cache.get('na_samples', []))
        calibrated = apply_calibration(img, geom, calibration)

    # PNG re-renders are the slow part of a cache-hit batch (each is
    # ~33 MB). Skip them if the output PNG already exists; pass
    # --rerender to force a refresh.
    norm_path = NORM_DIR / f'{stem}.png'
    annot_path = ANNOT_DIR / f'{stem}.png'
    NORM_DIR.mkdir(parents=True, exist_ok=True)
    if rerender or not norm_path.exists():
        save_png16(calibrated, norm_path)

    metrics = quantify_spots(img, geom, spots)

    ANNOT_DIR.mkdir(parents=True, exist_ok=True)
    if rerender or not annot_path.exists():
        overlay = annotate(calibrated, geom, spots, metrics, stem,
                           na_samples=na_samples)
        if not cv2.imwrite(str(annot_path), overlay):
            raise IOError(f'failed writing annotated PNG for {stem}')

    print(f'[done]   {stem}  '
          f'card=({geom.card_cx:.0f},{geom.card_cy:.0f}) r_card={geom.r_card:.0f}  '
          f'dish=({geom.dish_cx:.0f},{geom.dish_cy:.0f}) r_dish={geom.r_dish:.0f}')
    for m in metrics:
        print(f'         {m.sample:<10} sec={m.sector}  '
              f'r_eff={m.r_eff:5.1f}  area={m.area_px:6d}  '
              f'OD_chrom={m.OD_chrom:+.4f}')
    return metrics, meta


# ── Batch ──────────────────────────────────────────────────────

def _list_dngs() -> list:
    if not DATA_DIR.exists():
        sys.exit(f'data dir not found: {DATA_DIR}')
    return sorted(p for p in DATA_DIR.iterdir() if p.suffix.lower() == '.dng')


def main() -> None:
    cr_figure.apply_rcparams()

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--redetect', action='store_true',
                    help='Force re-detection (ignore cache)')
    ap.add_argument('--only', type=str, default=None,
                    help='Process a single image by stem (e.g. 48h_LB_1)')
    ap.add_argument('--plots-only', action='store_true',
                    help='Skip per-image processing; regenerate plots from CSVs')
    ap.add_argument('--panels', action='store_true',
                    help='Also write the standalone colony, levels, and effect panels')
    ap.add_argument('--skip-existing', action='store_true',
                    help='Skip any plate whose annotated PNG already exists; '
                         'reuse its rows from the prior quantification.csv. '
                         'Use after reshooting a small number of plates to '
                         'avoid re-decoding / re-annotating the rest.')
    ap.add_argument('--rerender', action='store_true',
                    help='Force re-rendering of normalized + annotated PNGs '
                         'even when they already exist.')
    ap.add_argument('--no-colony-strip', action='store_true',
                    help='Drop the colony crops from the effect-estimation '
                         'figure (overrides SHOW_COLONY_STRIP).')
    args = ap.parse_args()

    for d in (OUT_DIR, NORM_DIR, ANNOT_DIR, CACHE_DIR, FAILURES_DIR,
              PLOTS_DIR):
        d.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        rows = read_summary()
    else:
        dngs = _list_dngs()
        if args.only:
            dngs = [p for p in dngs if p.stem == args.only]
            if not dngs:
                sys.exit(f'no DNG matches stem {args.only!r}')
        if not dngs:
            sys.exit(f'no DNGs in {DATA_DIR}')

        rows = []
        failures = []
        single = args.only is not None
        existing_quant = load_existing_quant() if args.skip_existing else {}
        for path in dngs:
            if (args.skip_existing
                    and (ANNOT_DIR / f'{path.stem}.png').exists()
                    and path.stem in existing_quant):
                rows.extend(existing_quant[path.stem])
                print(f'[skip]   {path.stem} (annotated; reused '
                      f'{len(existing_quant[path.stem])} rows)')
                continue
            try:
                metrics, meta = process_image(path, redetect=args.redetect,
                                              rerender=args.rerender)
            except Exception as e:
                if single:
                    raise RuntimeError(f'while processing {path.name}: {e}') from e
                failures.append((path.name, str(e)))
                print(f'[FAIL]   {path.stem}: {e}')
                save_failure_diagnostic(path, str(e))
                continue
            # Success: clear any stale failure diagnostic from a prior run.
            prev_fail = FAILURES_DIR / f'{path.stem}.png'
            if prev_fail.exists():
                prev_fail.unlink()
            for m in metrics:
                rows.append({'image': path.stem, **meta, **asdict(m)})

        if failures:
            print(f'\n{len(failures)} plate(s) failed (excluded from CSV / plots):')
            for name, err in failures:
                print(f'  {name}: {err}')

        if single:
            return  # don't overwrite CSV / plots when running a single image
        if not rows:
            sys.exit('no plates succeeded; aborting')
        write_csvs(rows)

    plot_all(rows, PLOTS_DIR,
             show_strip=False if args.no_colony_strip else None,
             panels=args.panels)


if __name__ == '__main__':
    main()
