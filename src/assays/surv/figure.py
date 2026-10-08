"""
Survival curve figures for the SURV assay.

One killing-curve panel per treatment, plus a combined figure across all
three.

Curves are drawn in log10 space — the same space the ANOVA and the slope
models work in — so a plotted mean and error bar correspond directly to the
quantities being tested. Any timepoint containing a non-detect is drawn with
an open marker so an imputed value is never mistaken for a measurement.

These figures carry no title and no on-image caption: the treatment and its
exposure condition, the detection limit, the open-marker convention and the
error-bar definition all belong in the manuscript figure legend, not baked
into the export. In the combined panel the panel letter is what the legend
names each treatment by.
"""

from __future__ import annotations

import viz as fs  # selects the Agg backend on import
import matplotlib.pyplot as plt
import numpy as np

from . import config


def apply_rcparams():
    """Apply the house plotting style. Called by the driver."""
    fs.use()


# ── Aggregation ───────────────────────────────────────────────

def summarise_for_plot(df):
    """Mean and SEM of Log10_CFU per Treatment / Strain / Time.

    Aggregating in log10 space keeps the error bar symmetric in the space it
    is computed in, and matches the model the p-values come from.
    """
    grouped = (df.groupby(["Treatment", "Strain", "Time"], observed=True)
                 .agg(mean_log10=("Log10_CFU", "mean"),
                      sd_log10=("Log10_CFU", "std"),
                      n=("Log10_CFU", "size"),
                      n_censored=("Censored", "sum"))
                 .reset_index())
    grouped["sem_log10"] = grouped["sd_log10"] / np.sqrt(grouped["n"])
    return grouped


# ── Panels ────────────────────────────────────────────────────

def _draw_treatment(ax, summary, raw, treatment):
    """Draw one treatment's killing curves onto an existing axis."""
    panel = summary[summary["Treatment"] == treatment]
    panel_raw = raw[raw["Treatment"] == treatment]

    for strain in config.FIGURE_STRAIN_ORDER:
        series = panel[panel["Strain"].astype(str) == strain]
        if series.empty:
            continue
        colour = config.STRAIN_COLOURS.get(strain, fs.INK_SECONDARY)
        marker = config.STRAIN_MARKERS.get(strain, "o")
        ax.errorbar(series["Time"], series["mean_log10"],
                    yerr=series["sem_log10"], marker=marker,
                    color=colour, ecolor=colour,
                    elinewidth=fs.ERRORBAR_LW, markeredgecolor=colour,
                    label=config.STRAIN_LABELS.get(strain, strain), zorder=3)

        # Hollow marker where a timepoint contains a non-detect: the mean and
        # error bar there are partly or wholly set by the imputed value.
        censored = series[series["n_censored"] > 0]
        if not censored.empty:
            ax.plot(censored["Time"], censored["mean_log10"], linestyle="none",
                    marker=marker, markerfacecolor=fs.SURFACE,
                    markeredgecolor=colour, markeredgewidth=1.0, zorder=4)

    exposure = config.TREATMENT_EXPOSURE.get(treatment, treatment)
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel(f"Survival under {exposure}\n(log$_{{10}}$ CFU/mL)")
    ax.set_xticks(sorted(panel_raw["Time"].unique()))
    fs.finish_axis(ax)


def plot_survival_curves(df, output_dir=None, summary=None):
    """One killing-curve figure per treatment. Returns the paths written."""
    output_dir = output_dir or config.OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarise_for_plot(df) if summary is None else summary
    written = []

    for treatment in config.TREATMENT_ORDER:
        if treatment not in set(df["Treatment"]):
            continue
        fig, ax = plt.subplots(figsize=config.FIG_SIZE_SINGLE)
        _draw_treatment(ax, summary, df, treatment)
        fs.legend(ax, loc="upper right", title="Strain")
        fig.tight_layout()

        written.extend(fs.save(fig, output_dir / f"survival_{treatment.lower()}"))

    return written


def plot_survival_panel(df, output_dir=None, summary=None):
    """All treatments side by side, for the supplement.

    Both axes are shared: treatments run to different final timepoints, and
    letting each panel set its own x-range would rescale time between panels
    and distort the apparent death rate in a figure whose whole purpose is
    cross-treatment comparison.
    """
    output_dir = output_dir or config.OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = summarise_for_plot(df) if summary is None else summary

    treatments = [t for t in config.TREATMENT_ORDER if t in set(df["Treatment"])]
    fig, axes = plt.subplots(1, len(treatments),
                             figsize=config.FIG_SIZE_PANEL,
                             sharey=True, sharex=True)
    axes = np.atleast_1d(axes)

    for ax, treatment in zip(axes, treatments):
        _draw_treatment(ax, summary, df, treatment)

    # Shared ticks across every timepoint in the assay, so the panels share
    # one x range.
    ax0 = axes[0]
    ax0.set_xticks(sorted(df["Time"].unique()))
    ax0.set_xlim(df["Time"].min() - 3, df["Time"].max() + 3)

    # One legend for the figure, built from every strain drawn in any panel.
    # It rides in the last panel's own upper right: a killing curve descends
    # left to right, so that corner is empty in every panel, and EtOH stops
    # two timepoints short of the shared x range.
    handles, labels = fs.collect_legend(axes)
    fs.legend(axes[-1], handles, labels, loc="upper right", title="Strain")

    fs.tag_panels(axes)
    fig.tight_layout()
    return fs.save(fig, output_dir / "survival_all_treatments")[0]
