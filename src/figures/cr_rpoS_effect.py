"""Congo Red: colonies over strain levels, the rpoS effect drawn between them.

Two panels, both drawn by `assays.cr.panels`, so this figure and the assay's
standalone panels cannot disagree about how a level or an effect looks.

Each panel is one row of a single gridspec column, so the solver gives them one
left and right edge between them, and the strip draws its crops in the levels
panel's own data coordinates — which is what puts a crop over its column with
nothing moved after the fact.

The effect rides in the levels panel as a connector between the two columns it
is taken over. Its length is the paired estimate, on the shared absolute axis,
so the four effects are comparable to each other by eye; `assays.cr.figure`
keeps the per-dish differences and their intervals as a figure of their own.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import viz as fs  # selects the Agg backend on import

from assays.cr import panels
from assays.cr.config import RNG_SEED, SHOW_COLONY_STRIP


def plot(rows: list, comparisons: list, out_path: Path,
         times: list = None,
         show_strip: bool = None) -> None:
    """Colony strip over strain levels, the rpoS effect between the pairs.

    A shows each column's representative colony twice: as it photographs, and
    as the per-pixel `OD_chrom` field the metric is built on. B shows every
    strain's `OD_chrom` as its own distribution under an estimate rule, colour
    carrying the organism and fill the genotype, and carries the effect as a
    connector between each pair. Both stand over one set of columns, named
    once above A.

    Nothing joins a control to its mutant: they are separate isolates co-plated
    on one dish, so the dish is a block, not a repeated measure. The block
    still earns the statistic — control and knockout `OD_chrom` correlate
    across dishes at r ≈ 0.6 (every cell positive), so differencing within a
    dish tightens the estimate 1.3-1.6x at this reading over treating the
    strains as independent. The connector carries that paired Hodges-Lehmann
    estimate, with the Bonferroni-corrected signed-rank test of each effect
    against zero.

    `times` selects the readings to draw (default `TIMES`); with more than one
    they are grouped, since plate numbering is per reading and no plate is
    shared across that boundary. B scales to its data.
    """
    if show_strip is None:
        show_strip = SHOW_COLONY_STRIP
    conditions, cells, layout = panels.levels_context(
        rows, comparisons, times, key=show_strip)
    if not conditions or not comparisons:
        print('[skip] effect estimation: nothing to draw')
        return

    layers, scale = (panels.build_strip(comparisons, conditions, cells)
                     if show_strip else ([], None))

    fig, page = fs.stack(fs.COL_DOUBLE, [
        panels.strip_rows(layers),
        {'levels': 2.7},
    ])

    # One generator for the whole figure, so the jitter is reproducible.
    rng = np.random.default_rng(RNG_SEED)

    strip = [fig.add_subplot(page[f'crop{i}']) for i in range(len(layers))]
    for ax, layer in zip(strip, layers):
        panels.draw_crops(ax, layout, layer)
    if scale is not None:
        panels.draw_od_key(strip[-1], scale, layout)

    ax_levels = fig.add_subplot(page['levels'])
    panels.draw_levels(ax_levels, comparisons, cells, conditions, layout, rng)
    panels.draw_pair_effects(ax_levels, comparisons, cells, conditions, layout)

    # One key, over the topmost panel. Both panels share their columns and the
    # page is read downwards, so naming them here names them for both and
    # nothing below is met before it has a name. Only n stays with B, where the
    # marks it counts are.
    panels.label_columns(strip[0] if strip else ax_levels, comparisons,
                         conditions, layout, side='top')
    panels.label_counts(ax_levels, layout,
                        panels.column_counts(comparisons, conditions, cells))

    fs.fit(fig)
    # `box`, so the key over the strip is left unlettered: it names the
    # columns both panels stand on, not panel A's alone.
    fs.tag_panels(strip[:1] + [ax_levels], align=True, anchor='box')
    # After the letters, which freeze the layout: the crops need the axes'
    # final shape, and changing an extent moves nothing.
    panels.square_crops(strip)
    written = fs.save(fig, out_path)
    print('wrote ' + ', '.join(str(p) for p in written))
