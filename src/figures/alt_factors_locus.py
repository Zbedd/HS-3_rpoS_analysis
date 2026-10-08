"""HS-3 RpoD, RpoH, RpoE and RpoN, shown as loci and proteins.

Eight panels in four rows, one row per factor, each giving that factor the
two views the manuscript figure gives rpoS — panel B's neighbourhood across
the four genomes, then panel C's domain architecture:

    A  rpoD synteny        B  RpoD domains
    C  rpoH synteny        D  RpoH domains
    E  rpoE synteny        F  RpoE domains
    G  rpoN synteny        H  RpoN domains

Nothing is drawn here. The two panel functions come from
`seq.rpoS_locus.figure`, which the standalone panels also call, so a panel
printed on its own and the same panel on this page are one piece of code.

The panel includes the alternative factors selected in the manuscript configuration.
"""
from __future__ import annotations

from pathlib import Path

import viz as fs  # selects the Agg backend on import

from seq.rpoS_locus import domains, figure as locus_figure, io, panels

# One stacked row per σ factor, heights in inches. Each row carries a header
# and two lettered panels side by side; `stack` spaces the rows and reserves
# each letter's band, so nothing here adds padding of its own.
#
# The plot box is what the four headers cost: the page has to stay printable
# at its own size, since a supplement scaled down to fit loses the point of
# sizing type in journal column widths.
ROW_HEIGHT = 1.65

# The header's own band. It is a row of the same panel as the plots under it,
# so it sits one ROW_GAP above them and a full panel gap below the row before
# — which is what makes it read as their header rather than as a footer to the
# row above. Tall enough that the letters, which reach into the bottom of this
# band, stay clear of the text sitting at its top.
HEADER_HEIGHT = 0.30

# How a row divides. The synteny panel takes the wider share: it draws four
# genomes' worth of gene arrows over six kilobases and a column of row labels,
# where the domain cartoon draws four proteins and a key.
LOCUS_SPLIT = [0.57, 0.43]


def plot(out_path: Path, workdir: Path = None, arches: dict = None) -> list:
    """Draw the figure. Returns the paths written.

    `arches` is the InterProScan architecture map; passing the one the driver
    already loaded keeps the page from re-reading the cache.
    """
    cfg, paths = io.load(workdir)
    order = cfg["alt_factors"]["order"]
    if arches is None:
        arches = domains.ensure_scan(cfg, paths)

    names = cfg["alt_factors"]["sigma_names"]
    fig, page = fs.stack(
        fs.COL_DOUBLE,
        [{f"{key}_header": HEADER_HEIGHT, key: ROW_HEIGHT} for key in order],
        pad=0.20)

    axes = []
    for key in order:
        panels.draw_row_title(fig.add_subplot(page[f"{key}_header"]),
                              f"{fs.italic(key)} ({names[key]})")
        left, right = fs.split(page[key], LOCUS_SPLIT)
        ax_syn = fig.add_subplot(left)
        locus_figure.draw_alt_synteny(ax_syn, cfg, paths, key)
        ax_dom = fig.add_subplot(right)
        locus_figure.draw_alt_domains(ax_dom, cfg, paths, key, arches)
        axes += [ax_syn, ax_dom]

    fs.fit(fig)
    fs.tag_panels(axes, labels="ABCDEFGH", align=True)
    return fs.save(fig, out_path)
