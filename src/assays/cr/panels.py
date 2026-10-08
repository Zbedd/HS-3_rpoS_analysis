"""CR panel bodies — each draws into an axes the caller owns.

Nothing here creates a figure or writes a file, so the same panel serves the
assay's standalone figures and the composed figure in `figures`.

`level_positions` is the single owner of the levels x layout. The levels
panel, the colony strip above it and the effect panel below it all have to
agree on where a column sits, so the layout is computed once and handed to
each of them rather than derived three times. Every panel over those columns
draws in the same data coordinates, which is what lets a figure stack them
without moving anything afterwards.
"""
from __future__ import annotations

from typing import NamedTuple

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patheffects as mpatheffects
import matplotlib.transforms as mtransforms
import numpy as np

import multiplicity

import viz as fs  # selects the Agg backend on import

from .config import EPS, MEDIA_DISPLAY, PLOT_MEDIA, RNG_SEED, TIMES
from .io import cached_r_eff, load_colony_crops, load_existing_quant
from .stats import estimation_cells, signed_rank_p, wilcoxon_hl_ci


# Short axis label. The full definition of the metric is too long for an
# axis at column width, so it lives in the methods text — see the
# `quantify_spots` docstring for the formula.
OD_CHROM_LABEL = 'OD$_\\mathrm{chrom}$'


class EffectLayout(NamedTuple):
    """Where every effect marker sits, in the effect axis' data coordinates."""
    xpos: dict   # (comparison label, condition index) -> x
    xlim: tuple


def effect_positions(comparisons: list, n_conditions: int,
                     span: float = 0.17) -> EffectLayout:
    """Lay the effect series out: each dodged symmetrically about its
    condition tick, closer to each other than to the neighbouring conditions.

    Owned here so the markers, their intervals and any annotation over them
    agree on one x rather than deriving it three times.
    """
    k = len(comparisons)
    offsets = np.linspace(-span, span, k) if k > 1 else np.zeros(1)
    return EffectLayout({(comp['label'], i): i + dx
                         for comp, dx in zip(comparisons, offsets)
                         for i in range(n_conditions)},
                        (-0.5, n_conditions - 0.5))


def effect_positions_over(comparisons: list, conditions: list,
                          levels: 'LevelLayout') -> EffectLayout:
    """Effect positions that line up with the levels columns they summarise.

    The levels panel spaces its columns by strain and the effect panel by
    condition, so left alone the two disagree about where anything sits. They
    share a box, so this solves the effect axis' limits and dodge from the
    levels layout: a condition lands on its group centre, and a series on the
    midpoint of the two columns it is the difference of.
    """
    centres = [(lo + hi) / 2 for lo, hi in levels.group_span.values()]
    if len(centres) < 2:
        return effect_positions(comparisons, len(conditions))

    lo, hi = levels.xlim
    f0, f1 = ((c - lo) / (hi - lo) for c in (centres[0], centres[-1]))
    width = (len(centres) - 1) / (f1 - f0)

    # How far a comparison's own columns sit from their group centre, in the
    # levels axis' units, carried into the effect axis' units.
    media, time = conditions[0]
    offsets = [sum(levels.xpos[(comp['label'], media, time, n)]
                   for n in (comp['denominator'], comp['numerator'])) / 2
               - centres[0] for comp in comparisons]
    span = max(abs(o) for o in offsets) * width / (hi - lo)
    left = -f0 * width
    return EffectLayout(
        effect_positions(comparisons, len(conditions), span).xpos,
        (left, left + width))


# The estimate rule, shared by both panels: they show the same estimator of
# the same quantity, so they draw it the same way. It is neutral because hue
# is already spoken for — it names the organism, on the dishes the rule sits
# over — and a rule crossing a cloud of its own colour is what the mark it
# replaces could not do. The halo keeps it clear of them. `_` takes its length
# from `s` in points, so the rule prints the same width in both panels, whose
# x scales differ.
# 27 pt is 60% of the levels panel's 45 pt column pitch — the width the rule
# had when it was drawn in that panel's data units, before the effect panel
# (240 pt to the unit) had to share it.
RULE_SIZE = 729      # s is (length in points)², so 27 pt
RULE_LW = 1.0
RULE_HALO = [mpatheffects.withStroke(linewidth=2.4, foreground=fs.SURFACE)]

DOT_SIZE = 9         # one dish, in either panel

# Descending, so the first tier a p clears is the one it earns. n floors the
# reachable p at 2/2^n and the correction multiplies that by the number of
# comparisons, so the smallest value this figure can report is 4 x 2/2^12 =
# 0.002 and a '***' tier could never be earned — see `stats.signed_rank_p`.
SIG_TIERS = ((0.01, '**'), (0.05, '*'))


def stars_for(p: float, tiers: tuple = SIG_TIERS) -> str:
    """Star string for `p`; 'n.s.' when it clears no tier."""
    if not np.isfinite(p):
        return ''
    for cut, mark in tiers:
        if p < cut:
            return mark
    return 'n.s.'


def draw_effect(ax, comparisons: list, panels: list, labels: list,
                layout: EffectLayout, ylim: tuple = None,
                rng: np.random.Generator = None, label_x: bool = True) -> None:
    """Draw the rpoS effect onto `ax`: one dodged series per comparison, over
    the conditions named by `labels`.

    A rule is the Hodges-Lehmann estimate, the error bar the exact Wilcoxon
    signed-rank 95% CI (so 'CI clears the baseline' matches the signed-rank
    test), the dots the individual per-dish differences. The dots carry the
    organism's colour and the rule is neutral, exactly as in the levels panel:
    a difference belongs to an organism, and the estimate of it is the same
    statistic whichever organism that is.
    `panels` is one list of `estimation_cells` cells per comparison, aligned
    to `labels`. `ylim` fixes the value axis; without it the axis grows to
    hold the estimates, whiskers and dots. `label_x` names each series and each
    condition on the x axis; a figure that names the columns over a panel above
    this one passes False, since every marker sits on the pair it summarises.

    Shared by the standalone lollipop and the composed figure's effect panel,
    so the two cannot drift apart.
    """
    rng = rng if rng is not None else np.random.default_rng(RNG_SEED)
    baseline = 0.0
    ax.axhline(baseline, color=fs.AXIS, lw=0.7, zorder=1)

    all_dots, all_err, all_cen, legend_marks = [], [], [], []
    for comp, cells in zip(comparisons, panels):
        tone = _tone(comp)
        cen = np.asarray([c['hl'] for c in cells])
        lows = np.asarray([c['lo'] for c in cells])
        highs = np.asarray([c['hi'] for c in cells])
        # NaN CI (sample too small) -> zero-length whisker, marker still shown.
        err_lo = np.where(np.isfinite(lows), np.maximum(0, cen - lows), 0.0)
        err_hi = np.where(np.isfinite(highs), np.maximum(0, highs - cen), 0.0)
        gx = np.array([layout.xpos[(comp['label'], i)]
                       for i in range(len(labels))])

        for x, cell in zip(gx, cells):
            d = cell['diffs']
            if d.size == 0:
                continue
            ax.scatter(x + (rng.random(d.size) - 0.5) * 0.11, d,
                       s=DOT_SIZE, color=tone, alpha=0.9,
                       edgecolor='none', zorder=3)
            all_dots.append(d)
        ax.errorbar(gx, cen, yerr=[err_lo, err_hi], fmt='none',
                    ecolor=fs.INK, lw=fs.ERRORBAR_LW, zorder=4)
        ax.scatter(gx, cen, s=RULE_SIZE, marker='_', color=fs.INK,
                   linewidth=RULE_LW, zorder=5, path_effects=RULE_HALO)
        legend_marks.append(comp['label'])
        all_cen.append(cen)
        all_err.append(err_lo)
        all_err.append(err_hi)

    # Each series is named under its own marker and each condition under its
    # pair, the way the levels panel names its columns — a legend would repeat
    # what the axis can say in place.
    ax.set_xlim(*layout.xlim)
    if not label_x:
        ax.set_xticks([])
        fs.finish_axis(ax)
    else:
        order = sorted(layout.xpos, key=layout.xpos.get)
        ax.set_xticks([layout.xpos[k] for k in order])
        ax.set_xticklabels([k[0] for k in order])
        ax.tick_params(axis='x', pad=_TICK_PAD)
        fs.finish_axis(ax)
        # The same two tiers the levels columns carry, over the series each
        # condition holds rather than over columns.
        spans = [(min(xs) - 0.10, max(xs) + 0.10) for xs in
                 ([layout.xpos[k] for k in layout.xpos if k[1] == i]
                  for i in range(len(labels)))]
        _draw_spans(ax, spans, labels, _TICK_PAD + _LINE_PT + _SPAN_GAP,
                    'bottom', fs.SIZE_BODY, 'bold')

    # The effects are all negative, so the upper right is free for the key.
    if legend_marks:
        _series_legend(ax, comparisons)

    if ylim is not None:
        ax.set_ylim(*ylim)
        return
    finite = (np.concatenate([m[np.isfinite(m)] for m in all_cen])
              if all_cen else np.array([]))
    dots = np.concatenate(all_dots) if all_dots else np.array([])
    errs = np.concatenate(all_err) if all_err else np.array([])
    if finite.size:
        spread_cen = float(np.max(np.abs(finite - baseline))) * 1.4
        err_max = float(np.max(errs)) * 1.2 if errs.size else 0.0
        spread_dots = (float(np.max(np.abs(dots - baseline))) * 1.1
                       if dots.size else 0.0)
        span_y = max(0.15, spread_cen, err_max, spread_dots)
        ax.set_ylim(baseline - span_y, baseline + span_y)
    else:
        ax.set_ylim(baseline - 0.5, baseline + 0.5)


def draw_significance(ax, comparisons: list, panels: list,
                      layout: EffectLayout, tiers: tuple = SIG_TIERS) -> None:
    """Mark each effect against zero, above its own interval.

    The mark reads the corrected p: every cell the panel draws is one test in
    one family, so the correction runs over all of them together.
    """
    cells_drawn = [(comp, i, cell) for comp, cells in zip(comparisons, panels)
                   for i, cell in enumerate(cells) if cell['diffs'].size]
    adjusted = multiplicity.correct(
        [signed_rank_p(cell['diffs']) for _c, _i, cell in cells_drawn])
    for (comp, i, cell), p in zip(cells_drawn, adjusted):
        mark = stars_for(p, tiers)
        if mark in ('', 'n.s.'):
            continue
        # Above the dishes as well as the interval: several cells carry a
        # difference well above their own CI, and a mark placed on the
        # interval alone lands inside the cloud.
        top = max(cell['hi'] if np.isfinite(cell['hi']) else cell['hl'],
                  float(np.max(cell['diffs'])))
        ax.annotate(mark, xy=(layout.xpos[(comp['label'], i)], top),
                    xytext=(0, 4), textcoords='offset points',
                    ha='center', va='bottom', fontsize=fs.SIZE_SMALL,
                    color=fs.INK)


# Blank either side of the outermost column, in column units, and the extra
# right margin the OD-map key rides in — inside the axes, so a figure carrying
# the key keeps its printed width.
_PAD_X = 0.7
_KEY_SPAN = 0.45


class LevelLayout(NamedTuple):
    """Where every levels column sits, in the levels axis' data coordinates."""
    xpos: dict        # (comparison label, media, time, strain) -> x
    group_span: dict  # (media, time) -> (first x, last x)
    x_max: float      # the last column
    key: bool         # reserve the right margin for the OD-map key

    @property
    def xlim(self) -> tuple:
        """The x range every panel over these columns shares."""
        return -_PAD_X, self.x_max + _PAD_X + (_KEY_SPAN if self.key else 0.0)


# The three gaps, which must stay in this order: a pair is closer than two
# organisms, which are closer than two media, so the nesting reads before the
# key above it is. The pair gap also has to clear the jitter either side of it
# and leave the effect label room between the columns it belongs to — widening
# it alone does not help, since the pair gap is most of the axis and the points
# per unit fall with it. `_JITTER` is what buys that room.
_COL_STEP = 1.25     # a control to its knockout
_ORG_STEP = 0.45     # extra, one organism's pair to the next
_MEDIA_STEP = 1.00   # extra, one medium's block to the next
_JITTER = 0.30       # full width of the scatter around a column


def level_positions(comparisons: list, conditions: list,
                    key: bool = False) -> LevelLayout:
    """Lay the levels columns out.

    Columns run control, knockout for each background in turn, the
    backgrounds set slightly apart, the conditions further apart.
    """
    xpos, group_span = {}, {}
    x = 0.0
    for media, time in conditions:
        start = x
        for ci, comp in enumerate(comparisons):
            if ci:
                x += _ORG_STEP
            for name in (comp['denominator'], comp['numerator']):
                xpos[(comp['label'], media, time, name)] = x
                x += _COL_STEP
        group_span[(media, time)] = (start, x - _COL_STEP)
        x += _MEDIA_STEP
    return LevelLayout(xpos, group_span, max(xpos.values()), key)


def levels_context(rows: list, comparisons: list, times: list = None,
                   key: bool = False) -> tuple:
    """The conditions, their estimation cells and the column layout.

    The three things every panel over these columns needs, derived once so a
    standalone panel and a composed figure cannot disagree about any of them.
    """
    times = list(times) if times else list(TIMES)
    conditions = [(m, t) for t in times for m in PLOT_MEDIA]
    return (conditions, estimation_cells(rows, comparisons, conditions),
            level_positions(comparisons, conditions, key))


def condition_labels(conditions: list) -> list:
    """One label per condition, naming the reading only when several are drawn.

    Plate numbering is per reading, so two readings are two blocks of columns
    and the label has to say which.
    """
    single = len({t for _media, t in conditions}) == 1
    return [MEDIA_DISPLAY[m] if single else f'{MEDIA_DISPLAY[m]}\n{t}'
            for m, t in conditions]


def draw_levels(ax, comparisons: list, panels: list, conditions: list,
                layout: LevelLayout, rng: np.random.Generator) -> None:
    """Each strain's `OD_chrom` as its own jittered distribution under its
    estimate, control beside knockout.

    The estimate is drawn exactly as the effect panel draws its own — same
    rule, same width, same interval — because it is the same estimator of the
    same quantity, and a reader should not have to be told so. Colour carries
    the organism and fill the genotype, both on the dishes; an effect panel's
    dot keeps the colour and drops the fill, being a difference. A hairline
    joins the two spots that shared a dish — a statement about provenance,
    not a transformation: they are separate isolates co-plated on one plate,
    so the dish is a block, not a repeated measure. `label_columns` names the
    columns.
    """
    xpos = layout.xpos
    for comp, cells in zip(comparisons, panels):
        num, den = comp['numerator'], comp['denominator']
        for (media, time), cell in zip(conditions, cells):
            drawn = {}
            for name, vals in ((den, [d for _i, d, _n in cell['triples']]),
                               (num, [n for _i, _d, n in cell['triples']])):
                vals = np.asarray(vals, dtype=np.float64)
                key = (comp['label'], media, time, name)
                if not vals.size:
                    continue
                px = xpos[key]
                face, edge = fs.fill_for(name, fs.lineage_colour(name))
                est, lo, hi = wilcoxon_hl_ci(vals)
                jx = px + (rng.random(vals.size) - 0.5) * _JITTER
                drawn[name] = (jx, vals)
                ax.scatter(jx, vals, s=DOT_SIZE, facecolor=face,
                           edgecolor=edge, linewidth=0.7, alpha=0.9, zorder=3)
                if np.isfinite(lo) and np.isfinite(hi):
                    ax.errorbar(px, est, yerr=[[est - lo], [hi - est]],
                                fmt='none', ecolor=fs.INK,
                                lw=fs.ERRORBAR_LW, zorder=4)
                ax.scatter(px, est, s=RULE_SIZE, marker='_', color=fs.INK,
                           linewidth=RULE_LW, zorder=5,
                           path_effects=RULE_HALO)
            # `triples` is one row per dish, so index i is the same plate in
            # both columns.
            if den in drawn and num in drawn:
                (dx, dv), (nx, nv) = drawn[den], drawn[num]
                for seg in zip(zip(dx, dv), zip(nx, nv)):
                    ax.plot(*zip(*seg), color=fs.GRID, lw=0.5, zorder=2)

    ax.set_xlim(*layout.xlim)
    ax.set_ylabel(OD_CHROM_LABEL)
    fs.finish_axis(ax)
    _levels_legend(ax, comparisons)


# The pair connector. Its top sits on the control's estimate rule and its
# length is the paired estimate itself — not the gap between the two rules,
# which is a different statistic (HL is not linear, so HL of the differences
# is not the difference of the HLs; they disagree by up to 10% of the effect
# here). The length drawn is the number written beside it.
# An arrow, not a capped line: the intervals beside it are capped lines, and a
# drop is not an interval.
_LINK_LW = 0.9
# `annotate` scales an arrowstyle by the font size unless told otherwise, so
# the head is sized here in points and the scale pinned.
_LINK_HEAD = 5.0
_LINK_TEXT_PAD = 2.0  # from the arrow to its label, in points


def draw_pair_effects(ax, comparisons: list, panels: list, conditions: list,
                      layout: LevelLayout, tiers: tuple = SIG_TIERS) -> None:
    """The rpoS effect, drawn between the two columns it is taken over.

    One connector per pair, labelled with the fold change on the chromatic
    ratio. `OD_chrom` is a log ratio, so the difference the connector draws is
    already a log fold change and exp() carries it — and its interval — onto a
    multiplicative scale exactly, exp being monotone.

    Every pair the panel draws is one test in one family, so the marks read
    the corrected p, as the effect panel's do.
    """
    drawn = [(comp, media, time, cell)
             for comp, cells in zip(comparisons, panels)
             for (media, time), cell in zip(conditions, cells)
             if cell['diffs'].size]
    adjusted = multiplicity.correct(
        [signed_rank_p(cell['diffs']) for _c, _m, _t, cell in drawn])

    for (comp, media, time, cell), p in zip(drawn, adjusted):
        den, num = comp['denominator'], comp['numerator']
        control = np.asarray([d for _i, d, _n in cell['triples']], dtype=float)
        if not control.size:
            continue
        top = wilcoxon_hl_ci(control)[0]
        delta = float(cell['hl'])
        x = (layout.xpos[(comp['label'], media, time, den)]
             + layout.xpos[(comp['label'], media, time, num)]) / 2.0

        ax.annotate('', xy=(x, top + delta), xytext=(x, top),
                    arrowprops=dict(arrowstyle='-|>', color=fs.INK,
                                    linewidth=_LINK_LW,
                                    mutation_scale=_LINK_HEAD,
                                    shrinkA=0, shrinkB=0,
                                    path_effects=RULE_HALO),
                    zorder=6, annotation_clip=False)
        mark = stars_for(p, tiers)
        # Two lines, so the label fits the gap the columns leave it.
        ax.annotate(f'\u00d7{np.exp(delta):.2f}\n{mark}',
                    xy=(x, top + delta / 2.0),
                    xytext=(_LINK_TEXT_PAD, 0), textcoords='offset points',
                    ha='left', va='center', fontsize=fs.SIZE_TICK,
                    linespacing=1.15, color=fs.INK, zorder=7,
                    path_effects=RULE_HALO)


def _tone(comp: dict) -> str:
    """Colour of the organism a comparison is drawn for."""
    return comp.get('color') or fs.lineage_colour(comp['numerator'])


# Both panels key the same two channels, so both keys are built from the same
# handles: the levels panel carries colour and fill, the effect panel only
# colour, and the half they share cannot drift apart.
_KEY_MARK = 4.0


def _organism_handles(comparisons: list) -> list:
    """One dot per organism, in that organism's colour."""
    return [mlines.Line2D([], [], linestyle='none', marker='o',
                          markersize=_KEY_MARK, markeredgecolor='none',
                          markerfacecolor=_tone(comp), label=comp['label'])
            for comp in comparisons]


def _genotype_handles(comparisons: list) -> list:
    """The control open, the knockout solid, in neutral ink.

    Neutral because fill is the only thing these say; colour is the entries
    above them.
    """
    first = comparisons[0]
    return [
        mlines.Line2D([], [], linestyle='none', marker='o',
                      markersize=_KEY_MARK,
                      markerfacecolor=fs.SURFACE if name == first['denominator']
                      else fs.INK_SECONDARY,
                      markeredgecolor=fs.INK_SECONDARY, markeredgewidth=0.7,
                      label=fs.label(name))
        for name in (first['denominator'], first['numerator'])]


def _series_legend(ax, comparisons: list) -> None:
    """The effect panel's key: colour only, since a difference has no genotype.

    Keyed on the dots because they are what hue is on — the rule over them is
    the same statistic for every series and names none of them. It repeats the
    top half of the levels panel's key, so a panel read on its own still says
    what its colours mean.
    """
    if not comparisons:
        return
    fs.legend(ax, handles=_organism_handles(comparisons), loc='upper right')


def _levels_legend(ax, comparisons: list) -> None:
    """The levels panel's key: the organism on colour, the genotype on fill.

    Both channels are on the dishes here, so both belong in one key. The R2A
    block sits low, so the upper right is free.
    """
    if not comparisons:
        return
    fs.legend(ax, handles=_organism_handles(comparisons)
              + _genotype_handles(comparisons), loc='upper right')


def column_counts(comparisons: list, conditions: list, panels: list) -> dict:
    """How many dishes each column rests on, keyed like `LevelLayout.xpos`."""
    counts = {}
    for comp, cells in zip(comparisons, panels):
        for (media, time), cell in zip(conditions, cells):
            for name in (comp['denominator'], comp['numerator']):
                counts[(comp['label'], media, time, name)] = len(cell['diffs'])
    return counts


# What the key stacks beyond the axis it names, in points: the column names,
# then the organism span, then the condition span. The row need not reserve any
# of it — `viz.fit` gives the page back whatever the key takes.
_TICK_PAD = 2.5      # matplotlib's own, for a column name beside a spine
_HEADER_PAD = 8.0    # for one over an image, which has a hard edge to clear
_LINE_PT = 9.0
_SPAN_GAP = 4.5      # from the last line of text to the span rule beyond it
_SPAN_TICK = 3       # the rule's end ticks
_TIER_GAP = 4.0      # from one tier's label to the next tier's rule
_COLUMN_HALF = 0.30  # how far a span reaches either side of a column


def label_columns(ax, comparisons: list, conditions: list, layout: LevelLayout,
                  counts: dict = None, conditions_row: bool = True,
                  side: str = 'bottom') -> None:
    """Name every column on `ax`: the strain, its organism, its condition.

    The panel carrying the columns carries their names, so a figure stacking
    several such panels names them once, on whichever it likes. `side` is which
    edge they hang from — a stack is read downwards, so a key over the topmost
    panel reaches every panel under it before the reader meets any of them.

    Each tier spans what it names and is drawn the same way, so the nesting is
    legible before the words are: the condition over a whole block in a title's
    size and weight, the organism over each pair, the strain over its own
    column.

    `counts` adds each column's n to its label; `conditions_row` adds the
    condition spanning its whole block.
    """
    order = sorted(layout.xpos, key=layout.xpos.get)
    ax.set_xticks([])
    ax.set_xlim(*layout.xlim)

    tiers = [(
        [_span(layout.xpos[k]) for k in order],
        [fs.label(k[3]) + (f'\nn = {counts[k]}' if counts else '')
         for k in order],
        2 if counts else 1, fs.SIZE_TICK, 'normal',
    ), (
        [span for comp in comparisons
         for span in _pair_spans(comp, conditions, layout)],
        [comp['label'] for comp in comparisons for _ in conditions],
        1, fs.SIZE_TICK, 'normal',
    )]
    if conditions_row:
        tiers.append((
            [_span(lo, hi) for lo, hi in layout.group_span.values()],
            condition_labels(conditions), 1, fs.SIZE_BODY, 'bold',
        ))
    _stack_tiers(ax, tiers, _HEADER_PAD if side == 'top' else _TICK_PAD, side)


def _span(lo: float, hi: float = None) -> tuple:
    """The x range a tier's rule covers over one column, or over a run of
    them."""
    return lo - _COLUMN_HALF, (lo if hi is None else hi) + _COLUMN_HALF


def _stack_tiers(ax, tiers: list, pad: float, side: str) -> None:
    """Draw the key outwards from the axis, each tier clear of the one before.

    The stack is arithmetic on the type size, so a tier that gains a line
    pushes the rest out rather than landing on them.
    """
    drop = pad
    for spans, labels, lines, size, weight in tiers:
        _draw_spans(ax, spans, labels, drop, side, size, weight)
        drop += _TIER_GAP / 2 + _LINE_PT * lines + _SPAN_GAP


def label_counts(ax, layout: LevelLayout, counts: dict) -> None:
    """Each column's n, under the marks it counts."""
    order = sorted(layout.xpos, key=layout.xpos.get)
    ax.set_xticks([layout.xpos[k] for k in order])
    ax.set_xticklabels([f'n = {counts[k]}' for k in order])
    ax.set_xlim(*layout.xlim)


def _draw_spans(ax, spans: list, labels: list, drop: float, side: str,
                size: float, weight: str = 'normal') -> None:
    """A rule across each span of columns, with its name beyond it.

    One primitive for every tier of the key. Colour separates the four strains
    but does not say which two belong together; a span does, without spending a
    legend on it.
    """
    away = 1 if side == 'top' else -1
    edge = 1.0 if side == 'top' else 0.0
    trans = (mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
             + mtransforms.ScaledTranslation(0, away * drop / 72,
                                             ax.figure.dpi_scale_trans))
    for (lo, hi), text in zip(spans, labels):
        # End ticks as markers, so their length is in points and does not
        # follow the height of whichever panel carries the span.
        ax.plot([lo, hi], [edge, edge], transform=trans, color=fs.AXIS,
                lw=0.7, marker='|', markersize=_SPAN_TICK,
                markeredgewidth=0.7, clip_on=False, zorder=2)
        ax.annotate(text, xy=((lo + hi) / 2, edge), xycoords=trans,
                    xytext=(0, away * _TIER_GAP / 2),
                    textcoords='offset points', ha='center',
                    va='bottom' if side == 'top' else 'top',
                    fontsize=size, fontweight=weight, color=fs.INK,
                    annotation_clip=False)


def _pair_spans(comp: dict, conditions: list, layout: LevelLayout) -> list:
    """The x range each of a comparison's pairs of columns covers."""
    spans = []
    for media, time in conditions:
        xs = [layout.xpos[(comp['label'], media, time, n)]
              for n in (comp['denominator'], comp['numerator'])]
        spans.append((min(xs) - 0.30, max(xs) + 0.30))
    return spans


# Crop half-width as a multiple of the median colony radius, so a colony
# sits clear of its frame with some agar around it.
_CROP_PAD = 2.2


def _fetch_crops(comparisons: list, conditions: list, panels: list) -> tuple:
    """Raw linear crops for every cell's representative plate, keyed
    (comparison label, media, time, sample), with the plate each came from.

    One decode per plate; both comparisons can land on the same one.
    """
    radii, wanted = [], {}
    for comp, cells in zip(comparisons, panels):
        for (media, time), cell in zip(conditions, cells):
            if cell['rep'] is None:
                continue
            stem = f'{time}_{media}_{cell["rep"]}'
            samples = (comp['denominator'], comp['numerator'])
            wanted[(comp['label'], media, time)] = (stem, samples)
            radii.extend(r for r in (cached_r_eff(stem, s) for s in samples)
                         if r)
    if not radii:
        return {}, {}

    half_px = int(round(_CROP_PAD * float(np.median(radii))))
    by_stem = {}
    for stem, samples in wanted.values():
        by_stem.setdefault(stem, set()).update(samples)
    loaded = {stem: load_colony_crops(stem, sorted(names), half_px)
              for stem, names in by_stem.items()}

    crops, source = {}, {}
    for (label, media, time), (stem, samples) in wanted.items():
        for name in samples:
            img = loaded.get(stem, {}).get(name)
            if img is not None:
                crops[(label, media, time, name)] = img
                source[(label, media, time, name)] = (stem, name)
    return crops, source


def load_strip(comparisons: list, conditions: list, panels: list) -> dict:
    """Colony crops as they photograph, keyed like `LevelLayout.xpos`.

    One gamma and one stretch across the whole set, so crops stay comparable
    between columns.
    """
    crops, _ = _fetch_crops(comparisons, conditions, panels)
    if not crops:
        return {}
    crops = {k: np.clip(v, 0, 1) ** (1 / 2.2) for k, v in crops.items()}
    # Raw is linear light and sits in the bottom few percent of the range, so
    # display needs a black and white point, not just a scale.
    flat = np.concatenate([c.ravel() for c in crops.values()])
    lo, hi = (float(v) for v in np.percentile(flat, [2.0, 99.8]))
    if hi - lo > 1e-6:
        crops = {k: np.clip((v - lo) / (hi - lo), 0, 1)
                 for k, v in crops.items()}
    return crops


def load_od_maps(comparisons: list, conditions: list, panels: list) -> tuple:
    """Per-pixel `OD_chrom` over the same crops, and the norm to draw them by.

    The colony summary is a log of ratios of means; this is the mean-free
    per-pixel form, log(ref_G/G) − log(ref_R/R), so the two differ slightly by
    Jensen's inequality and only this one is displayed. Each spot is
    referenced to its own annulus, whose medians are already columns in
    `quantification.csv`, so nothing is re-detected.

    Agar sits at zero by construction, so the norm is centred and symmetric
    about it: a quarter of crop pixels fall just below zero as agar noise, and
    a one-sided norm would drive them to full saturation. The end comes from a
    robust percentile over every crop at once.
    """
    crops, source = _fetch_crops(comparisons, conditions, panels)
    if not crops:
        return {}, None
    quant = load_existing_quant()
    refs = {(stem, r['sample']): (r['ref_R'], r['ref_G'])
            for stem, rws in quant.items() for r in rws
            if isinstance(r.get('ref_R'), float)}

    maps = {}
    for key, img in crops.items():
        ref = refs.get(source[key])
        if ref is None:
            continue
        ref_R, ref_G = ref
        maps[key] = (np.log(ref_G / np.maximum(img[..., 1], EPS))
                     - np.log(ref_R / np.maximum(img[..., 0], EPS)))
    if not maps:
        return {}, None
    flat = np.concatenate([m.ravel() for m in maps.values()])
    v = max(float(np.percentile(flat, 99.5)), 1e-6)
    return maps, mcolors.Normalize(vmin=-v, vmax=v)


class StripLayer(NamedTuple):
    """One row of the colony strip: images keyed like `LevelLayout.xpos`, what
    the row shows, and how `imshow` should draw them."""
    images: dict
    label: str = ''
    imshow_kw: dict = {}


def build_strip(comparisons: list, conditions: list, panels: list) -> tuple:
    """The strip's layers, and the scale its OD-map layer is drawn by.

    The OD-map layer is drawn last, so `draw_od_key` belongs on the last axes.
    Returns ([], None) when neither the crops nor the reference medians are
    available.
    """
    layers, scale = [], None
    photo = load_strip(comparisons, conditions, panels)
    if photo:
        layers.append(StripLayer(photo, 'colony'))
    od, norm = load_od_maps(comparisons, conditions, panels)
    if od:
        cmap = fs.colormap('od_chrom')
        scale = mcm.ScalarMappable(norm=norm, cmap=cmap)
        layers.append(StripLayer(od, OD_CHROM_LABEL,
                                 {'cmap': cmap, 'norm': norm}))
    return layers, scale


# Crop half-width in column units. Columns sit 1.0 apart, so this is as wide
# as a crop goes without touching its neighbour. `CROP_ROW_IN` is the row
# height at which a crop that wide comes out square at double-column width.
CROP_HALF = 0.56
CROP_ROW_IN = 0.60


def strip_rows(layers: list) -> dict:
    """Row heights for a colony strip."""
    return {f'crop{i}': CROP_ROW_IN for i in range(len(layers))}


def draw_crops(ax, layout: LevelLayout, layer: StripLayer) -> None:
    """One strip layer: each column's crop, at that column's x.

    The images carry the levels panel's own coordinates, so a crop sits over
    its column with nothing moved afterwards and `label_columns` names them
    exactly as it names the levels panel's. `square_crops` gives them square
    pixels once the axes has its final size.
    """
    for key, x in layout.xpos.items():
        img = layer.images.get(key)
        if img is not None:
            ax.imshow(img, extent=(x - CROP_HALF, x + CROP_HALF, 0.0, 1.0),
                      aspect='auto', interpolation='nearest', **layer.imshow_kw)
    ax.set_xlim(*layout.xlim)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(sorted(layout.xpos.values()))
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.grid(False)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    if layer.label:
        ax.set_ylabel(layer.label, rotation=0, ha='right', va='center',
                      labelpad=5, fontsize=fs.SIZE_SMALL,
                      color=fs.INK_SECONDARY)


def square_crops(axes) -> None:
    """Give the crops square pixels, now that the axes have their final size.

    An image keeps its shape only in an axes free to change its own, which a
    panel sharing a box with the panels under it is not — so the extent is what
    gives. Each crop takes as much of its row as it can without touching its
    neighbour, and loses height rather than overlap.
    """
    for ax in axes:
        box = ax.get_window_extent()
        span = abs(np.diff(ax.get_xlim())[0])
        square = 0.5 * span * box.height / box.width
        half = min(CROP_HALF, square)
        height = half / square
        for image in ax.images:
            x0, x1, _lo, _hi = image.get_extent()
            mid = (x0 + x1) / 2
            image.set_extent((mid - half, mid + half,
                              0.5 - height / 2, 0.5 + height / 2))


def draw_od_key(ax, scale, layout: LevelLayout) -> None:
    """The OD-map colour scale, in the right margin the layout reserves for it.

    An inset tracks its parent, so the key follows the strip without the
    figure placing it.
    """
    lo, hi = layout.xlim
    left = (layout.x_max + CROP_HALF + 0.14 - lo) / (hi - lo)
    bar = ax.figure.colorbar(
        scale, cax=ax.inset_axes([left, 0.0, 0.10 / (hi - lo), 1.0]))
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=2, width=0.6, labelsize=fs.SIZE_SMALL,
                       colors=fs.INK_SECONDARY)
