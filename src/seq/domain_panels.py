"""Shared protein-domain panel drawing. Creates no figure and writes no file.

`figure.render_domain_cartoon` wraps `draw_domain_cartoon` for RQ2's own
supplementary figures; `figures.rpoS_locus_and_mutant` calls it directly for
the manuscript figure. One function, so the two cannot disagree about what a
domain architecture looks like.

Type sizes are a parameter because the two callers print at different widths:
RQ2's own figure is 11.5 in wide and keeps the sizes it was drawn at, while
the manuscript panel is a column of a 7.09 in page and takes the house ladder.
"""
from __future__ import annotations

import math

import matplotlib.patches as patches

import viz as fs

from .protein_domains import PF_R11, PF_R12, PF_R2, PF_R3, PF_R4, PF_R4_2

DOMAIN_COLOURS = fs.DOMAIN_COLOURS
UNMAPPED_COLOUR = fs.DOMAIN_UNMAPPED

DOMAIN_LABELS = {
    PF_R11:  "σ70 r1.1",
    PF_R12:  "σ70 r1.2",
    PF_R2:   "σ70 r2",
    PF_R3:   "σ70 r3",
    PF_R4:   "σ70 r4",
    PF_R4_2: "σ70 r4.2",
    "PF04546": "σ70 NCR",       # non-essential region; Group 1 factors only
    # σ54, which is not a σ70-family protein: a different fold with its own
    # activator-interacting, core-binding and DNA-binding domains.
    "PF00309": "σ54 AID",
    "PF04963": "σ54 core",
    "PF04552": "σ54 DBD",
}

# The sizes RQ2's supplementary figures are drawn at.
RQ2_TYPE = {"block": 7, "label": 9, "label_hs3": 9.5, "annotation": 7,
            "title": 11, "axis": None}

# The house ladder, for the manuscript page.
HOUSE_TYPE = {"block": fs.SIZE_SMALL, "label": fs.SIZE_SMALL,
              "label_hs3": fs.SIZE_TICK, "annotation": fs.SIZE_SMALL,
              "title": fs.SIZE_BODY, "axis": fs.SIZE_TICK}

HIGHLIGHT = "#fff8d6"     # the HS-3 rows
BACKBONE = "#eaeaea"      # the part of a protein no signature covers


def draw_domain_cartoon(ax, arches: list, *, title: str = None,
                        right_label_family: str = None,
                        use_full_pfams: bool = False,
                        type_scale: dict = None,
                        emphasis: str = "band",
                        label_mode: str = "block",
                        align_on: str = None,
                        shape: str = "block",
                        show_axis: bool = True,
                        scale_bar: tuple = None,
                        legend_ncol: int = None) -> float:
    """One row per protein, coloured blocks per Pfam location.

    `arches` is [(display label, DomainArchitecture, is_hs3)]. HS-3 rows take
    a bold label, a heavier backbone and a highlight band, so the candidate
    reads against its references at a glance. Returns the longest protein
    length, which sets the row width the caller may want for alignment.

    `right_label_family` selects the family reported in the right margin:
    "NCBIfam" surfaces the gene-specific curated HMMs (TIGR02394 rpoS_proteo,
    TIGR02393 RpoD_Cterm) that carry the σ⁷⁰ call; "Pfam" surfaces accessions,
    for a selected protein family. None omits it.

    `use_full_pfams` draws every Pfam-A signature rather than the selected family alone, for the unfiltered companion view.

    `label_mode` is where a domain is named. "block" writes the name inside
    every block, picking ink by contrast against that block's fill — which
    means the ink changes colour from block to block. "top" writes each name
    once, above the first row, in a single ink. "legend" moves the names out
    of the rows entirely into a key. All three are the secondary encoding the
    six-hue palette needs; the last two are one typeface colour rather than
    four.

    `show_axis` draws the residue axis. Turn it off and pass `scale_bar` as
    (span, label) for a cartoon that carries a scale instead of coordinates —
    which is what a panel aligned on a domain wants, since its axis measures
    distance from that domain rather than position in the protein.

    `align_on` is a Pfam accession every row is shifted to start at, instead
    of every row starting at residue 1. Aligning puts the domains in columns,
    so the comparison the panel exists to make is read down the page rather
    than inferred across ragged rows; the cost is that the axis becomes
    residues relative to that domain, not absolute position. A protein lacking
    the domain is left unshifted.

    `shape` draws each domain as a "block" or as an "arrow" — the same
    pentagon a gene takes on a chromosome panel, pointing N to C.

    `legend_ncol` wraps the key. It defaults to one row, which is right for a
    panel as wide as the page; a panel sharing a row with another needs the
    names over two or three lines instead, and the space that takes is added
    below the rows rather than taken out of them.

    `emphasis` is how an HS-3 row is set apart: "band" gives it the highlight
    strip and heavier outlines RQ2's supplementary figures use, "label" leaves
    every row drawn alike and carries the distinction on the bold row name
    alone — which is what a manuscript panel wants, where the highlight reads
    as a claim about the data rather than as a pointer.
    """
    sizes = dict(RQ2_TYPE if type_scale is None else type_scale)
    banded = emphasis == "band"
    named_above = label_mode == "top"
    keyed = label_mode == "legend"
    top_row = len(arches) - 1        # rows are drawn bottom-up
    seen = {}                        # accession -> (label, leftmost start)
    n = len(arches)
    max_len = max(a.length for _, a, _ in arches) if arches else 100

    # Ordered families take the whole ramp across what this panel draws: index
    # against the canonical six and a four-domain panel gets only the middle
    # of the range, where neighbouring steps are hardest to tell apart.
    present = {acc for _, arch, _ in arches
               for acc in (arch.pfams_full if use_full_pfams else arch.pfams)}
    ordered = [acc for acc in fs.SIGMA70_ORDER if acc in present]

    def block_colour(acc):
        if acc in ordered:
            return fs.ordinal(ordered.index(acc), len(ordered))
        return fs.domain_colour(acc)

    rows = list(reversed(arches))
    shifts = []
    for _, arch, _ in rows:
        spans = (arch.pfams.get(align_on) or arch.pfams_full.get(align_on)
                 if align_on else None)
        shifts.append(spans[0][0] if spans else 0)
    left_edge = min((-d for d in shifts), default=0)
    right_edge = max(arch.length - d for (_, arch, _), d in zip(rows, shifts))

    for i, (label, arch, is_hs3) in enumerate(rows):
        y = i
        dx = shifts[i]
        heavy = is_hs3 and banded
        if is_hs3 and banded:
            ax.add_patch(patches.Rectangle(
                (left_edge - max_len * 0.25, y - 0.42),
                (right_edge - left_edge) * 1.45, 0.84,
                facecolor=HIGHLIGHT, edgecolor="none", zorder=0))
        ax.add_patch(patches.Rectangle(
            (-dx, y - 0.18), arch.length, 0.36, facecolor=BACKBONE,
            edgecolor="black", linewidth=0.9 if heavy else 0.6, zorder=1))
        pfam_hits = arch.pfams_full if use_full_pfams else arch.pfams
        for acc, spans in pfam_hits.items():
            colour = block_colour(acc)
            lab = DOMAIN_LABELS.get(acc, acc)

            # One block per InterProScan location — no merging across gaps.
            for s, e in spans:
                s, e = s - dx, e - dx
                # The key follows the domains N to C, which is the order they
                # are read in, not the order HMMER happened to report them.
                if acc not in seen or s < seen[acc][1]:
                    seen[acc] = (lab, s)
                if shape == "arrow":
                    fs.arrow(ax, s, e, 1, y, 0.60, colour, edgecolor="black",
                             lw=0.9 if heavy else 0.6, zorder=2,
                             min_head=max_len * 0.035)
                else:
                    ax.add_patch(patches.Rectangle(
                        (s, y - 0.30), e - s, 0.60, facecolor=colour,
                        edgecolor="black", linewidth=0.9 if heavy else 0.6,
                        zorder=2))
                if keyed:
                    pass
                elif named_above:
                    if i == top_row:
                        ax.text((s + e) / 2, y + 0.42, lab, ha="center",
                                va="bottom", fontsize=sizes["block"],
                                color=fs.INK, weight="bold", zorder=3)
                else:
                    # Ink is picked by contrast against the block it sits on.
                    ax.text((s + e) / 2, y, lab, ha="center", va="center",
                            fontsize=sizes["block"], color=fs.ink_on(colour),
                            weight="bold", zorder=3)
        if right_label_family == "Pfam":
            right_label = ", ".join(sorted(arch.pfams.keys())) or "—"
        elif right_label_family == "NCBIfam":
            right_label = ", ".join(sorted(arch.nfams.keys())) or "—"
        else:
            right_label = None
        if right_label is not None:
            ax.text(right_edge + max_len * 0.02, y,
                    f"  {right_label_family}: {right_label}",
                    ha="left", va="center", fontsize=sizes["annotation"],
                    zorder=3)
        ax.text(left_edge - max_len * 0.02, y, label, ha="right",
                va="center",
                fontsize=sizes["label_hs3"] if heavy else sizes["label"],
                weight="bold" if is_hs3 else "normal", zorder=3)

    key_cols = min(legend_ncol or len(seen) or 1, len(seen) or 1)
    key_rows = math.ceil(len(seen) / key_cols) if keyed and seen else 1
    right_pad = 0.45 if right_label_family in ("NCBIfam", "Pfam") else 0.05
    ax.set_xlim(left_edge - max_len * 0.25, right_edge + max_len * right_pad)
    ax.set_ylim(-0.6 - (0.95 + 0.40 * (key_rows - 1)
                        if scale_bar or keyed else 0.0),
                n - 0.4 + (0.35 if named_above else 0.0))
    if show_axis:
        axis_label = ("residue position" if not align_on else
                      f"residues from {DOMAIN_LABELS.get(align_on, align_on)}")
        ax.set_xlabel(axis_label, **({} if sizes["axis"] is None
                                     else {"fontsize": sizes["axis"]}))
    else:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    if scale_bar:
        span, label = scale_bar
        fs.scale_bar(ax, left_edge, -0.95, span, label, tick=0.10)
    if keyed and seen:
        order = sorted(seen, key=lambda acc: seen[acc][1])
        handles = [patches.Polygon([(0, 0)], closed=True,
                                   facecolor=block_colour(acc),
                                   edgecolor="black", linewidth=0.6)
                   for acc in order]
        fs.legend(ax, handles, [seen[acc][0] for acc in order],
                  loc="lower right", ncol=key_cols, handlelength=1.1,
                  handletextpad=0.45, columnspacing=1.1, borderpad=0.0,
                  fontsize=sizes["block"])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=sizes["title"], loc="left", weight="bold")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    return max_len
