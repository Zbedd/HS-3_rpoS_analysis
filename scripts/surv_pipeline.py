"""
SURV driver — stress-survival killing curves for HS-3.

Loads data/SURV/survival_data.csv, runs the statistics, writes the figures,
and emits the tables and run manifest to outputs/SURV/.

    python scripts/surv_pipeline.py            # full run
    python scripts/surv_pipeline.py --help
"""

from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import scipy
import statsmodels

from assays.surv import config, io, stats
from assays.surv import figure  # selects the Agg backend on import


# ── Output helpers ────────────────────────────────────────────

def _write_table(table, path):
    """Write a DataFrame to CSV, reporting the row count."""
    table.to_csv(path, index=False, encoding="utf-8")
    print(f"  wrote {path.name:<40} {len(table):>4} rows")
    return path


def emit_manifest(df, results, censoring, output_dir, *, data_file,
                  lod_divisor, figures, sensitivity, panels=False):
    """Record the inputs, conventions and package versions behind this run."""
    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "package_versions": {
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "input_file": str(data_file),
        "output_dir": str(output_dir),
        "options": {
            "lod_divisor": lod_divisor,
            "figures": figures,
            "standalone_panels": panels,
            "sensitivity_sweep": sensitivity,
        },
        "n_observations": int(len(df)),
        "treatments": [t for t in config.TREATMENT_ORDER
                       if t in set(df["Treatment"])],
        "strains": [s for s in config.STRAIN_ORDER
                    if s in set(df["Strain"].astype(str))],
        "conventions": {
            "cfu_per_ml": "count * 10**Dilution_Factor / (Volume_plated_uL / 1000)",
            "lod": "10**Dilution_Factor / (Volume_plated_uL / 1000)",
            "imputed_value": "LOD / lod_divisor, for non-detects only",
            "model_space": "log10(CFU/mL)",
            "anova": "Log10_CFU ~ C(Strain) * C(Time), type II",
            "within_timepoint": "Tukey HSD fitted separately at each timepoint",
            "unadjusted_contrast": "Welch t-test between strains at a timepoint",
            "slopes": "linregress(Time, Log10_CFU) per strain, and per replicate",
            "slope_comparison": "pooled ANCOVA C(Strain):Time, and Welch t-test "
                                "on per-replicate slopes",
            "multiplicity": {
                "methods": config.CORRECTION_METHODS,
                "applied_to": "unadjusted p-values across all reported contrasts",
                "scope": "the repository corrects by Bonferroni throughout",
                "note": "Tukey p-values are already family-wise corrected "
                        "within their timepoint and are not re-adjusted.",
            },
            "alpha": config.ALPHA,
        },
        "censoring": {
            "n_censored": int(df["Censored"].sum()),
            "detail": censoring.to_dict(orient="records"),
            "excluded_fully_censored_groups":
                results["excluded_groups"].to_dict(orient="records"),
        },
        "replicate_diagnostic": results["replicates"].to_dict(orient="records"),
    }

    path = output_dir / "run_manifest.json"
    path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"  wrote {path.name}")
    return path


# ── CLI ───────────────────────────────────────────────────────

def _positive_float(value):
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Stress-survival (killing curve) analysis for HS-3.")
    parser.add_argument("--data", type=Path, default=None,
                        help="input CSV (default: data/SURV/survival_data.csv)")
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="output directory (default: outputs/SURV/)")
    parser.add_argument("--lod-divisor", type=_positive_float,
                        default=config.LOD_DIVISOR,
                        help="where to place non-detects below the detection "
                             f"floor (default: {config.LOD_DIVISOR}, i.e. "
                             f"LOD/{config.LOD_DIVISOR})")
    parser.add_argument("--no-figures", action="store_true",
                        help="skip figure rendering")
    parser.add_argument("--panels", action="store_true",
                        help="also write a standalone killing curve for each treatment")
    parser.add_argument("--no-sensitivity", action="store_true",
                        help="skip the LOD sensitivity sweep")
    args = parser.parse_args(argv)

    data_file = args.data or config.DATA_FILE
    output_dir = args.output_dir or config.OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/6] Loading survival counts")
    df = io.load_survival(path=args.data, lod_divisor=args.lod_divisor)
    censoring = io.censoring_summary(df)
    print(f"  {len(df)} plates | {int(df['Censored'].sum())} non-detects "
          f"| LOD/{args.lod_divisor:g} imputation")

    print("[2/6] Statistics")
    results = stats.run_stats(df)
    for name, table in results.items():
        print(f"  {name:<20} {len(table):>4} rows")
    if not results["excluded_groups"].empty:
        for _, row in results["excluded_groups"].iterrows():
            print(f"  excluded from contrasts: {row.Treatment} t={row.Time} "
                  f"(all {row.n_plates} plates are non-detects)")

    print("[3/6] Writing tables")
    _write_table(df, output_dir / "survival_processed.csv")
    _write_table(censoring, output_dir / "censoring_summary.csv")
    for name, filename in [
        ("anova", "anova.csv"),
        ("tukey", "within_timepoint_contrasts.csv"),
        ("slopes", "death_rate_slopes.csv"),
        ("replicate_slopes", "death_rate_slopes_by_replicate.csv"),
        ("slope_comparisons", "slope_comparisons.csv"),
        ("replicates", "replicate_diagnostic.csv"),
        ("excluded_groups", "excluded_groups.csv"),
    ]:
        _write_table(results[name], output_dir / filename)

    if args.no_sensitivity:
        print("[4/6] LOD sensitivity sweep — skipped")
    else:
        print("[4/6] LOD sensitivity sweep")
        sweep = stats.lod_sensitivity(
            lambda lod_divisor: io.load_survival(path=args.data,
                                                 lod_divisor=lod_divisor),
            primary=args.lod_divisor)
        for name, table in sweep.items():
            _write_table(table, output_dir / f"sensitivity_{name}.csv")

    if args.no_figures:
        print("[5/6] Figures — skipped")
    else:
        print("[5/6] Figures")
        figure.apply_rcparams()
        summary = figure.summarise_for_plot(df)
        if args.panels:
            for path in figure.plot_survival_curves(df, output_dir, summary=summary):
                print(f"  wrote {path.name}")
        print(f"  wrote {figure.plot_survival_panel(df, output_dir, summary=summary).name}")

    print("[6/6] Manifest")
    emit_manifest(df, results, censoring, output_dir,
                  data_file=data_file, lod_divisor=args.lod_divisor,
                  figures=not args.no_figures,
                  sensitivity=not args.no_sensitivity,
                  panels=args.panels and not args.no_figures)
    print(f"\nDone. Outputs in {output_dir}")


if __name__ == "__main__":
    main()
