"""The CR panels as standalone figures — one panel each.

The manuscript figure lives in `figures.cr_rpoS_effect` and composes the
colony strip and levels panel, with effect arrows, so a panel cannot say one thing alone and another in
company. What does differ is only what a panel must carry when it has no
neighbours: standing alone it names its own conditions, and the effect panel
spaces its series evenly rather than under levels columns that are not there.

Figures draw from `PLOT_SAMPLES` and `PLOT_MEDIA`, never `SAMPLES` /
`MEDIA_ORDER` — the latter include conditions kept only for CSV provenance.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import viz as fs  # selects the Agg backend on import

from . import panels
from .config import RNG_SEED, RPOS_MINUS


def apply_rcparams() -> None:
    """Apply the house plotting style. Idempotent — call from `main()`."""
    fs.use()


def _save(fig, out_path: Path) -> None:
    written = fs.save(fig, out_path)
    print('wrote ' + ', '.join(str(p) for p in written))


def plot_colonies(rows: list, comparisons: list, out_path: Path,
                  times: list = None) -> None:
    """Each column's representative colony, as it photographs and as the
    per-pixel `OD_chrom` field the metric is built on."""
    conditions, cells, layout = panels.levels_context(
        rows, comparisons, times, key=True)
    layers, scale = panels.build_strip(comparisons, conditions, cells)
    if not layers:
        print('[skip] colonies: no crops available')
        return
    fig, page = fs.stack(fs.COL_DOUBLE, [panels.strip_rows(layers)])
    axes = [fig.add_subplot(page[f'crop{i}']) for i in range(len(layers))]
    for ax, layer in zip(axes, layers):
        panels.draw_crops(ax, layout, layer)
    panels.label_columns(axes[-1], comparisons, conditions, layout)
    if scale is not None:
        panels.draw_od_key(axes[-1], scale, layout)
    fs.fit(fig)
    panels.square_crops(axes)
    _save(fig, out_path)


def plot_levels(rows: list, comparisons: list, out_path: Path,
                times: list = None) -> None:
    """Every strain's `OD_chrom` as its own distribution, control beside
    knockout, under a Hodges-Lehmann estimate and its exact Wilcoxon 95% CI."""
    conditions, cells, layout = panels.levels_context(rows, comparisons, times)
    fig, page = fs.stack(fs.COL_DOUBLE, [{'levels': 2.0}])
    ax = fig.add_subplot(page['levels'])
    panels.draw_levels(ax, comparisons, cells, conditions, layout,
                       np.random.default_rng(RNG_SEED))
    panels.label_columns(ax, comparisons, conditions, layout,
                         counts=panels.column_counts(comparisons, conditions,
                                                     cells))
    fs.fit(fig)
    _save(fig, out_path)


def plot_effect(rows: list, comparisons: list, out_path: Path,
                ylabel: str = None, ylim: tuple = None,
                times: list = None) -> None:
    """The per-dish rpoS effect: Hodges-Lehmann estimate, exact Wilcoxon
    signed-rank 95% CI, and a mark for the signed-rank test against zero."""
    conditions, cells, layout = panels.levels_context(rows, comparisons, times)
    fig, page = fs.stack(fs.COL_DOUBLE, [{'effect': 2.3}])
    ax = fig.add_subplot(page['effect'])
    effect = panels.effect_positions(comparisons, len(conditions))
    panels.draw_effect(ax, comparisons, cells,
                       panels.condition_labels(conditions), effect,
                       ylim=ylim, rng=np.random.default_rng(RNG_SEED))
    panels.draw_significance(ax, comparisons, cells, effect)
    ax.set_ylabel(ylabel if ylabel is not None else
                  f'{panels.OD_CHROM_LABEL}  ({RPOS_MINUS} − WT)')
    fs.fit(fig)
    _save(fig, out_path)
