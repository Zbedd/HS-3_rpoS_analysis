# SURV — migration map

Where the `HS3_survival` code landed. Source repository:
<https://github.com/icalimlim3/HS3_survival> (upstream HEAD `e669766`).
All 9 upstream commits are preserved in this repository's history.

| Upstream | Here |
|---|---|
| `survival_stats.py` | `src/assays/surv/stats.py` |
| `survival_plots.py` | `src/assays/surv/figure.py` |
| `Survival_plots/survival_plotter.py` | folded into `src/assays/surv/figure.py` |
| `survival_data.xlsx` | `data/SURV/survival_data.xlsx` (+ `.csv`) |
| `Survival_plots/output*.png` | regenerated into `outputs/SURV/` |

The manuscript workflow starts from the recorded plate counts in
`data/SURV/survival_data.csv`, with the source workbook retained alongside it.
See [the data notes](../../../data/SURV/README.md) for the input schema and
analysis conventions.

## Function map

| Upstream | Here |
|---|---|
| `run_stats(file_path)` | `stats.run_stats(df)` — takes the loaded frame; ANOVA split into `run_anova`, Tukey into `run_tukey` |
| `calculate_death_slopes(df, treatment_name)` | `stats.calculate_death_slopes` — returns a frame with slope CIs; the pairwise slope test moved to `stats.compare_death_slopes`, and `stats.replicate_slopes` adds the per-culture fit |
| inline CFU/LOD arithmetic | `io.calculate_cfu_ml`, `io.calculate_lod`, `io.load_survival` |
| inline plotting block | `figure.plot_survival_curves`, `figure.plot_survival_panel` |

## Layout

```
assays.surv.config    paths, strain order, LOD convention, plot style
assays.surv.io        load_survival(): CFU/mL, LOD, censoring flag
assays.surv.stats     ANOVA, within-timepoint contrasts, slopes, multiplicity
assays.surv.figure    killing curves
scripts/surv_pipeline.py driver — python scripts/surv_pipeline.py --help
```

`io.load_survival` is the only definition of CFU/mL, the detection limit and
the censoring flag; both the statistics and the figures consume the frame it
returns. It refuses input it cannot analyse cleanly — unknown strain or
treatment labels, non-numeric dilutions or volumes — rather than letting the
affected rows disappear from some tables and not others.

## Outputs

`outputs/SURV/` (untracked):

| File | Contents |
|---|---|
| `survival_processed.csv` | every plate with CFU/mL, LOD, imputed value, censoring flag |
| `anova.csv` | two-way ANOVA per treatment, plus Levene and an interaction flag |
| `within_timepoint_contrasts.csv` | strain contrasts at each timepoint, three p-value columns |
| `death_rate_slopes.csv` | slope per strain with 95% CI |
| `death_rate_slopes_by_replicate.csv` | slope per biological replicate |
| `slope_comparisons.csv` | pooled ANCOVA and replicate-level slope tests |
| `replicate_diagnostic.csv` | variance attributable to replicate identity |
| `censoring_summary.csv`, `excluded_groups.csv` | non-detects, and cells dropped for being wholly censored |
| `sensitivity_*.csv` | every family repeated across the LOD conventions |
| `survival_all_treatments.{png,pdf}` | combined supplementary panel (default) |
| `survival_{heat,h2o2,etoh}.{png,pdf}` | individual treatments, with `--panels` |
| `run_manifest.json` | inputs, options, conventions, package versions |
