"""How panels are arranged. `style` owns how they look.

A composed figure is a stack of lettered panels, each one or more named rows
tall:

    fig, page = fs.stack(fs.COL_DOUBLE, [{'levels': 2.2}, {'effect': 3.0}])
    ax_levels = fig.add_subplot(page['levels'])
    ax_effect = fig.add_subplot(page['effect'])
    fs.tag_panels([ax_levels, ax_effect], align=True)

Rows are sized in inches, so one keeps its printed height when a neighbour
changes, which ratios alone cannot promise. The gaps between them are not
given: they come from `style`. A row's height covers whatever its panel hangs
under the axis — tick labels, brackets, a condition row — because the solver
takes that space out of the axes rather than out of the row.

Every figure built here uses constrained layout. Nested gridspecs are the
normal case for a composed figure, and `tight_layout` does not solve them.
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from .style import INK, LETTER_BAND, PANEL_GAP, ROW_GAP, SIZE_TITLE

__all__ = ["figure", "stack", "fit", "split", "tag_panels", "collect_legend"]

# Two panels count as one column when their left edges land this close, in
# figure widths. Wider than solver jitter, narrower than any real gutter.
_COLUMN_TOL = 0.01


def figure(width: float, heights: dict, *, pad: float = 0.0,
           hspace: float = 0.30, **gridspec_kw):
    """A constrained-layout figure whose rows are `heights` inches tall.

    `heights` maps a row name to its height in inches; the figure is their
    sum plus `pad`, and the same numbers become the gridspec's ratios, so the
    two cannot drift apart. Returns the figure and a name -> `SubplotSpec`
    mapping. `stack` is the same thing with the gaps supplied for you.
    """
    if not heights:
        raise ValueError("figure() needs at least one row")
    names = list(heights)
    inches = [float(heights[n]) for n in names]
    fig = plt.figure(figsize=(width, sum(inches) + pad), layout="constrained")
    grid = fig.add_gridspec(len(inches), 1, height_ratios=inches,
                            hspace=hspace, **gridspec_kw)
    return fig, {name: grid[i] for i, name in enumerate(names)}


def stack(width: float, panels, *, pad: float = 0.0, **gridspec_kw):
    """A page of lettered panels, spaced from the style guide.

    `panels` is one dict per lettered panel, mapping each of its row names to
    that row's height in inches. Rows of one panel are separated by `ROW_GAP`,
    panels by `PANEL_GAP` and the band the next letter needs; the figure is
    their sum plus `pad`. Returns the figure and a name -> `SubplotSpec` map.

    A height is the panel's plot box, not its cell: whatever it hangs under
    its axis — tick labels, brackets, a condition row — is added by `fit`,
    which every stacked figure calls once its panels are drawn.

    Every row lands in one gridspec column, so the solver gives every panel the
    same left and right edge. A panel drawn in another's data coordinates then
    sits over the columns it labels with nothing moved afterwards.
    """
    panels = [p for p in panels if p]
    if not panels:
        raise ValueError("stack() needs at least one row")
    rows = {}
    for i, panel in enumerate(panels):
        for j, (name, height) in enumerate(panel.items()):
            if j:
                rows[f" gap {len(rows)}"] = ROW_GAP
            rows[name] = float(height)
        if i < len(panels) - 1:
            rows[f" gap {len(rows)}"] = PANEL_GAP + LETTER_BAND
    fig, cells = figure(width, rows, pad=pad, hspace=0.0, **gridspec_kw)
    fig._stack_rows = list(rows.values())
    return fig, {name: cell for name, cell in cells.items()
                 if not name.startswith(" ")}


def fit(fig, tries: int = 3) -> None:
    """Grow a stacked figure until each row's plot box is the height it asked.

    Constrained layout shares a gridspec's height ratios out over the axes, not
    over the cells holding them, so every decoration on the page shrinks every
    row by one common factor. Measuring that factor and giving the shortfall
    back to the page is what makes the inches `stack` was handed real ones —
    and what lets a panel declare its plot box and say nothing about the labels
    hanging under it. Call it once the panels are drawn.
    """
    rows = getattr(fig, "_stack_rows", None)
    if not rows:
        return
    for _ in range(tries):
        fig.canvas.draw()
        scale = None
        for ax in fig.axes:
            spec = ax.get_subplotspec()
            if spec is None or spec.get_gridspec().nrows != len(rows):
                continue
            scale = (ax.get_position().height * fig.get_figheight()
                     / rows[spec.rowspan.start])
            break
        if scale is None:
            return
        short = sum(rows) * (1.0 - scale)
        if short < 0.01:
            return
        fig.set_figheight(fig.get_figheight() + short)


def split(cell, into, *, axis: str = "x", **gridspec_kw) -> list:
    """Subdivide one cell. Returns the pieces in order.

    `into` is a count for equal pieces, or a sequence of ratios. `axis` is
    "x" for columns (the default) or "y" for rows.
    """
    count = into if isinstance(into, int) else len(into)
    ratios = None if isinstance(into, int) else list(into)
    if axis == "x":
        sub = cell.subgridspec(1, count, width_ratios=ratios, **gridspec_kw)
        return [sub[0, i] for i in range(count)]
    if axis == "y":
        sub = cell.subgridspec(count, 1, height_ratios=ratios, **gridspec_kw)
        return [sub[i, 0] for i in range(count)]
    raise ValueError(f"axis must be 'x' or 'y', not {axis!r}")


def tag_panels(axes, labels: str = "ABCDEFGH", x: float = -0.03,
               y: float = 1.02, size: float = None, align: bool = False,
               anchor: str = "tight") -> None:
    """Bold panel letters (A, B, C…) above the top-left of each axis.

    By default each letter is placed in its own axis's coordinates, which
    suits a row of panels of equal geometry. `align` instead measures the
    laid-out panels and shares one left edge down each column, so a panel
    carrying a wider y-axis label than the panel above it does not push its
    letter out of line. `x` and `y` are then offsets in inches from that edge
    and from each panel's top, and the layout is frozen: the letters are
    written at absolute figure coordinates, so nothing may move a panel after
    this. Call it last.

    `anchor` is what the letter sits above: "tight" clears whatever the panel
    draws over its axes, which is right for a title; "box" is the plot box
    itself, which leaves a key drawn over the topmost panel outside the
    lettering — a key shared by every panel belongs to none of them.
    """
    size = size or SIZE_TITLE
    axes = list(axes)
    if not align:
        for ax, tag in zip(axes, labels):
            ax.text(x, y, tag, transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=size, fontweight="bold", color=INK, clip_on=False)
        return

    fig = axes[0].figure
    fig.canvas.draw()
    fig.set_layout_engine("none")
    to_fig = fig.transFigure.inverted()
    # Tight extents include the y-axis label and tick labels, so a shared edge
    # clears the widest of them rather than the widest axes box.
    extents = [to_fig.transform(ax.get_tightbbox(fig.canvas.get_renderer()))
               for ax in axes]
    starts = [ax.get_position().x0 for ax in axes]
    tops = [ax.get_position().y1 if anchor == "box" else extent[1][1]
            for ax, extent in zip(axes, extents)]

    # Only panels in the same column share an edge. Panels side by side in one
    # row each keep their own, or their letters would stack on top of one
    # another. Columns are read off the laid-out axes boxes, not the gridspec,
    # so an axes moved by hand is still grouped by where it ended up.
    edges = {}
    for start, extent in zip(starts, extents):
        column = min((c for c in edges if abs(c - start) <= _COLUMN_TOL),
                     key=lambda c: abs(c - start), default=start)
        edges[column] = min(edges.get(column, extent[0][0]), extent[0][0])

    dx = x / fig.get_figwidth()
    dy = (y - 1.0) / fig.get_figheight()
    for tag, start, top in zip(labels, starts, tops):
        column = min(edges, key=lambda c: abs(c - start))
        fig.text(edges[column] + dx, top + dy, tag, ha="left", va="bottom",
                 fontsize=size, fontweight="bold", color=INK)


def collect_legend(axes) -> tuple:
    """One (handles, labels) pair covering every series drawn in any panel.

    First occurrence wins, so the order follows the panels as drawn.
    """
    handles, labels = [], []
    for ax in axes:
        for handle, text in zip(*ax.get_legend_handles_labels()):
            if text not in labels:
                handles.append(handle)
                labels.append(text)
    return handles, labels
