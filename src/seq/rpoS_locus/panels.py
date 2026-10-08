"""Panels for the rpoS-locus figure. Each draws into an axes the caller owns.

No function here creates a figure, writes a file, or sets an rcParam —
`viz` owns all of that, and the caller owns the page.

Gene arrows are drawn in base pairs and every panel is to scale, including
the 6 kb vector in panel D: an insertion that big is the point of the panel,
so it is not compressed behind an axis break.
"""
from __future__ import annotations

import math

import matplotlib.patches as patches

import viz as fs

# What a locus panel draws its subject in: the σ factor the panel is about.
# `viz` anchors the domain ramp on this same value, so a subject's gene arrow
# and the first block of its domain cartoon are one colour and the two panels
# identify the same RpoS protein in the locus and domain views.
#
# The supplement's four factors each take it in their own row: a row is about
# one factor, and no page puts two of them, or one of them and rpoS, together.
SUBJECT_FILL = fs.colour("rpoS")

# The conserved cassette. One hue per gene, and the same hue wherever that
# gene is drawn — panel A's mark, panel B's arrow, panel D's two copies and
# the standalone construct map's arm.
CASSETTE_ORDER = ("surE", "pcm", "nlpD", "rpoS")
CASSETTE_COLOURS = {
    "surE": fs.HUES["blue"],
    "pcm": fs.HUES["violet"],
    "nlpD": fs.HUES["yellow"],
    "rpoS": SUBJECT_FILL,
}

# pRE118 by what each feature does, so the map reads as four functions rather
# than seven arbitrary blocks, and so no vector hue repeats a cassette hue
# anywhere in the figure. sacB and its promoter are one unit; oriT and the tra
# genes are the conjugation machinery.
VECTOR_CLASS = {
    "KanR": "selection",
    "sacB": "counter-selection",
    "sacB promoter": "counter-selection",
    "R6K γ ori": "replication",
    "oriT": "transfer",
    "traJ": "transfer",
    "traK": "transfer",
}
# Selection takes the kanR strain's own colour: the cassette this class marks
# is what names that strain, so the two are one identity rather than a clash
# between the vector map and the strain key.
VECTOR_CLASS_COLOURS = {
    "selection": fs.colour("kanR"),
    "counter-selection": fs.HUES["green"],
    "replication": fs.HUES["red"],
    "transfer": fs.HUES["magenta"],
}

CONTEXT_FILL = fs.role("backdrop")   # drawn, but not part of the claim
VECTOR_FILL = fs.role("feature")     # a vector feature with no class of its own
GRID_RING = fs.role("ring")          # unannotated span of a circular replicon


# What a σ factor locus other than rpoS is drawn in. The factor itself takes
# `SUBJECT_FILL`, exactly as rpoS does in panel B, and its neighbours take the
# cassette fills in the order they are drawn. Read off the cassette rather than
# restated, so the supplement cannot drift from panel B.
NEIGHBOUR_FILLS = tuple(CASSETTE_COLOURS[gene] for gene in CASSETTE_ORDER
                        if gene != "rpoS")


def cassette_fill(key: str) -> str:
    """Colour for a cassette gene, the same in every panel that draws it."""
    return CASSETTE_COLOURS[key]


def locus_palette(order, anchor: str) -> dict:
    """One colour per core gene of a σ factor locus, keyed on ortholog group.

    Reserving `SUBJECT_FILL` for the anchor rather than for a position in the
    order is what lets the four supplement loci read alike: rpoD closes its
    operon, rpoE opens its own, and rpoN sits in the middle of one.
    """
    neighbours = [g for g in order if g != anchor]
    if len(neighbours) > len(NEIGHBOUR_FILLS):
        raise ValueError(f"{anchor}: {len(neighbours)} neighbours, "
                         f"{len(NEIGHBOUR_FILLS)} fills")
    fills = dict(zip(neighbours, NEIGHBOUR_FILLS))
    fills[anchor] = SUBJECT_FILL
    return fills


def vector_palette(labels) -> dict:
    """One colour per vector feature, keyed on what the feature does.

    Keyed on the label rather than on iteration order, so the integration
    and standalone construct maps agree and adding a feature cannot shuffle the others.
    """
    return {label: VECTOR_CLASS_COLOURS.get(VECTOR_CLASS.get(label),
                                            VECTOR_FILL)
            for label in dict.fromkeys(labels)}


def construct_palette(vector_features) -> dict:
    """Vector fills plus the arm, which takes the value rpoS has everywhere
    else — one gene, one colour, across four panels."""
    palette = vector_palette(f.label for f in vector_features)
    palette["rpoS'"] = cassette_fill("rpoS")
    return palette


def gene_symbol(protein_name: str) -> str:
    """`RpoS` -> `rpoS`. Only the first letter changes case; `RpoD` -> `rpoD`."""
    return protein_name[0].lower() + protein_name[1:]


def _prime(symbol: str) -> str:
    """A gene symbol carrying the prime that marks a partial copy."""
    return fs.italic(symbol) + "′"


def construct_name(vector_name: str) -> str:
    """Name of the vector carrying the cloned homology arm."""
    return f"{vector_name}::" + fs.italic("rpoS") + "′"


_arrow = fs.arrow


def _block(ax, start, end, y, height, facecolor, *, lw=0.5, edgecolor=None,
           zorder=2, **kwargs):
    ax.add_patch(patches.Rectangle(
        (start, y - height / 2), end - start, height, facecolor=facecolor,
        edgecolor=edgecolor or fs.INK, linewidth=lw, zorder=zorder, **kwargs))


def _bare(ax):
    """Strip every axis decoration. The cartoon panels carry their own scale
    bar or ring instead of a matplotlib axis."""
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    ax.patch.set_visible(False)
    ax.grid(False)


def _break_mark(ax, x0, x1, y, height):
    """The double squiggle that stands for sequence removed from the drawing.

    Everything either side of it is to scale; the distance across it is not.
    """
    span, steps = x1 - x0, 41
    amplitude = span * 0.13
    for centre in (x0 + span * 0.34, x0 + span * 0.66):
        points = [(centre + amplitude * math.sin(4 * math.pi * i / (steps - 1)),
                   y + (i / (steps - 1) - 0.5) * height * 1.7)
                  for i in range(steps)]
        ax.plot([x for x, _ in points], [v for _, v in points],
                color=fs.INK_SECONDARY, lw=0.9, zorder=6,
                solid_capstyle="round")


def _cassette_legend(ax, tracks, fills: dict) -> None:
    """Key for the conserved cassette, in the order it is drawn."""
    order = [g.ortholog for g in tracks[0].genes if g.ortholog]
    handles = [patches.Polygon([(0, 0)], closed=True,
                               facecolor=fills[key], edgecolor=fs.INK,
                               linewidth=0.5) for key in order]
    fs.legend(ax, handles, [fs.italic(key) for key in order],
              loc="lower right", ncol=len(order), handlelength=1.1,
              handletextpad=0.45, columnspacing=1.1, borderpad=0.0,
              fontsize=fs.SIZE_SMALL)


_scale_bar = fs.scale_bar


def _radial_label(ax, angle, radius, text, *, inward: bool = False, **kwargs):
    """Text at a point on a ring, set clear of it.

    A label outside the ring reads away from the centre; one inside reads back
    towards it, so `inward` flips which side of the anchor the text extends to.
    """
    x, y = radius * math.cos(angle), radius * math.sin(angle)
    outer = "left" if x > 0.02 else ("right" if x < -0.02 else "center")
    flipped = {"left": "right", "right": "left", "center": "center"}
    ha = flipped[outer] if inward else outer
    va = "center" if abs(x) > 0.02 else ("bottom" if y > 0 else "top")
    ax.text(x, y, text, ha=ha, va=va, **kwargs)


# ── Panel A: the chromosome and its σ factors ─────────────────

def draw_chromosome(ax, length: int, sigmas: list, *, focus: str = "RpoS",
                    housekeeping: str = "RpoD", note: str = None) -> None:
    """The closed chromosome as a ring, each σ factor marked at its locus.

    `sigmas` is [(protein name, Feature)]. Position 0 is at the top and
    coordinates run clockwise. The focus and housekeeping factors take the
    rpoS and WT strain hues; the rest are neutral, so the panel reads as
    the two highlighted factors among the five shown.
    """
    _bare(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-1.58, 1.58)
    ax.set_ylim(-1.26, 1.26)
    radius, width = 1.0, 0.070

    def angle(bp):
        return math.pi / 2 - 2 * math.pi * bp / length

    ax.add_patch(patches.Wedge((0, 0), radius, 0, 360, width=width,
                               facecolor=GRID_RING, edgecolor=fs.AXIS,
                               linewidth=0.5, zorder=1))

    ax.text(0, 0, f"{length / 1e6:.2f} Mb", ha="center", va="center",
            fontsize=fs.SIZE_TICK, color=fs.INK_SECONDARY)

    # rpoS takes the hue it carries throughout the figure, so the ring ties to
    # the panels below it. Every other sigma factor is neutral: naming them is
    # the point, and multiple hues here would collide with the cassette.
    highlight = {focus: cassette_fill("rpoS"), housekeeping: fs.INK_SECONDARY}
    for name, feature in sigmas:
        a = angle(feature.midpoint)
        colour = highlight.get(name, fs.INK_SECONDARY)
        lead = name in highlight
        # A CDS is far too small to see at this scale, so the mark is a fixed
        # angular width rather than the gene drawn to scale.
        centre = math.degrees(a)
        span = 2.2 if lead else 1.4
        ax.add_patch(patches.Wedge((0, 0), radius, centre - span, centre + span,
                                   width=width, facecolor=colour,
                                   edgecolor="none", zorder=3))
        ax.plot([radius * math.cos(a), (radius + 0.085) * math.cos(a)],
                [radius * math.sin(a), (radius + 0.085) * math.sin(a)],
                color=colour, lw=1.4 if lead else 0.9, zorder=3)
        _radial_label(ax, a, radius + 0.13, fs.italic(gene_symbol(name)),
                      fontsize=fs.SIZE_TICK if lead else fs.SIZE_SMALL,
                      fontweight="bold" if lead else "normal",
                      color=fs.INK if lead else fs.INK_SECONDARY)
    if note:
        ax.text(0, -1.22, note, ha="center", va="bottom",
                fontsize=fs.SIZE_SMALL, color=fs.INK_MUTED)


# ── A row header, over a pair of panels ───────────────────────

def draw_row_title(ax, text: str) -> None:
    """The σ factor a row of the supplement is about, over both its panels.

    Set at the top of the axes it is given rather than centred in it: the
    panel letters below reach up into the same band, and this keeps the two
    from meeting.
    """
    _bare(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 1.0, text, ha="center", va="top", fontsize=fs.SIZE_BODY,
            fontweight="bold", color=fs.INK)


# ── Panel B: synteny across the comparators ───────────────────

def draw_synteny(ax, tracks: list, *, flanks: tuple, anchor: str = "rpoS",
                 fills: dict = None, name_locals: bool = False,
                 ribbons: bool = True) -> None:
    """One row per genome, every window anchored and oriented on `anchor`.

    Core genes are coloured and joined between adjacent rows; context genes
    are drawn but not labelled, since nothing here tested them.

    `name_locals` writes the symbol a genome gives a core gene where it
    differs from the group key — algU over P. aeruginosa's rpoE. Off for the
    manuscript panel, where the one such case (K-12 calls surE umpG) is not
    what that panel is about.
    """
    _bare(ax)
    fills = CASSETTE_COLOURS if fills is None else fills
    height, gap = 0.36, 1.15
    ys = {t.tag: -i * gap for i, t in enumerate(tracks)}
    anchor_len = max(g.end for g in tracks[0].genes if g.ortholog == anchor)

    lo, hi = -flanks[0], anchor_len + flanks[1]
    # The margin on the left is the row labels' column: genes are clamped to
    # the window so a clipped arrow cannot run underneath a label.
    ax.set_xlim(lo - 0.115 * (hi - lo), hi)
    ax.set_ylim(ys[tracks[-1].tag] - 1.30, gap * 0.30)

    if ribbons:
        for upper, lower in zip(tracks, tracks[1:]):
            bottom = {g.ortholog: g for g in lower.genes if g.ortholog}
            # Walked in track order, not set order: overlapping translucent
            # ribbons paint in the sequence they are drawn, and a set would
            # reorder them from run to run.
            for a in upper.genes:
                key = a.ortholog
                if not key or key not in bottom:
                    continue
                b = bottom[key]
                ax.add_patch(patches.Polygon(
                    [(a.start, ys[upper.tag] - height / 2),
                     (a.end, ys[upper.tag] - height / 2),
                     (b.end, ys[lower.tag] + height / 2),
                     (b.start, ys[lower.tag] + height / 2)],
                    closed=True, facecolor=fills[key], alpha=0.30,
                    edgecolor="none", zorder=0))

    for track in tracks:
        y = ys[track.tag]
        ax.plot([lo, hi], [y, y], color=fs.AXIS, lw=0.6, zorder=1)
        for gene in track.genes:
            start, end = max(gene.start, lo), min(gene.end, hi)
            if end - start < 1:
                continue
            colour = fills[gene.ortholog] if gene.ortholog else CONTEXT_FILL
            _arrow(ax, start, end, gene.strand, y, height, colour,
                   edgecolor=fs.INK if gene.ortholog else fs.INK_MUTED,
                   lw=0.5 if gene.ortholog else 0.4, min_head=260)
        # Context genes are named on every row: they are what differs between
        # the genomes, so leaving them anonymous would leave the reader taking
        # the conserved order on trust.
        for gene in track.genes:
            # A gene clipped by the window keeps its arrow but loses its name:
            # a label floating past the edge of the row names nothing.
            if not lo < gene.midpoint < hi:
                continue
            local = gene.local if name_locals else ""
            name = (local or gene.ortholog) if gene.ortholog else gene.label
            if gene.pseudo and name:
                # ψ is the standard mark, and it is set upright: the symbol
                # is not part of the gene name.
                text, ink = "ψ" + fs.italic(name), fs.INK
            elif gene.ortholog:
                if not local or local == gene.ortholog:
                    continue
                text, ink = fs.italic(local), fs.INK
            elif gene.label:
                text, ink = fs.italic(gene.label), fs.INK_MUTED
            else:
                continue
            ax.text(gene.midpoint, y - height / 2 - 0.09, text,
                    ha="center", va="top", fontsize=fs.SIZE_SMALL, color=ink)
        ax.text(lo - 0.014 * (hi - lo), y, track.label, ha="right", va="center",
                fontsize=fs.SIZE_TICK, color=fs.INK,
                style="normal" if track.label == "HS-3" else "italic",
                fontweight="bold" if track.label == "HS-3" else "normal")

    _scale_bar(ax, lo + 0.02 * (hi - lo), ys[tracks[-1].tag] - 0.70,
               1000, "1 kb")
    _cassette_legend(ax, tracks, fills)


# ── Panel D: the loop-in, designed and as sequenced ───────────

_TIER_Y = {"construct": 2.30, "wt": 1.15, "mutant": 0.0}

def draw_integration(ax, layout: dict, *, palette: dict = None,
                     vector_name: str = "pRE118") -> None:
    """Design over outcome: the construct, the target, and what resulted.

    All three tiers share one base-pair scale, so the 6 kb the vector adds is
    read off the page rather than taken on trust.
    """
    _bare(ax)
    height = 0.34
    arm_lo, arm_hi = layout["arm"]
    palette = palette or {}

    # The construct tier is not collapsed, so it is the widest thing on the
    # page and sets the right edge.
    lo = min(s.start for tier in ("construct", "wt", "mutant")
             for s in layout[tier])
    hi = max(s.end for tier in ("construct", "mutant") for s in layout[tier])
    pad = 0.03 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-0.92, 3.16)

    def draw_tier(name, segments, above=False, gap=None):
        y = _TIER_Y[name]
        sign = 1 if above else -1
        va = "bottom" if above else "top"
        edges = [min(s.start for s in segments), max(s.end for s in segments)]
        spans = ([(edges[0], gap[0]), (gap[1], edges[1])] if gap
                 else [(edges[0], edges[1])])
        for x0, x1 in spans:
            ax.plot([x0, x1], [y, y], color=fs.AXIS, lw=0.6, zorder=1)
        for seg in segments:
            if seg.kind in ("vector", "backbone"):
                _block(ax, seg.start, seg.end, y, height, VECTOR_FILL,
                       edgecolor=fs.INK_SECONDARY)
                if seg.kind == "backbone":
                    ax.text((seg.start + seg.end) / 2, y, vector_name,
                            ha="center", va="center", fontsize=fs.SIZE_SMALL,
                            color=fs.INK_SECONDARY, zorder=4)
            elif seg.kind == "arm":
                _block(ax, seg.start, seg.end, y, height,
                       cassette_fill("rpoS"), edgecolor=fs.INK)
            elif seg.kind == "marker":
                _block(ax, seg.start, seg.end, y, height,
                       palette.get(seg.label, VECTOR_FILL), edgecolor=fs.INK)
            elif seg.kind == "truncated":
                # A square end, not an arrow head: the head would assert a
                # stop codon this copy does not have.
                _block(ax, seg.start, seg.end, y, height,
                       cassette_fill("rpoS"), edgecolor=fs.INK)
            else:
                _arrow(ax, seg.start, seg.end, seg.strand, y, height,
                       cassette_fill(seg.label)
                       if seg.label in CASSETTE_ORDER else CONTEXT_FILL,
                       edgecolor=fs.INK if seg.label in CASSETTE_ORDER
                       else fs.INK_MUTED, min_head=240)
            if seg.kind in ("gene", "truncated", "arm", "marker"):
                # Vector-borne genes keep the label the plasmid map gives them,
                # so panels D and E name the same block the same way.
                text = (_prime(seg.label.rstrip("′'"))
                        if seg.kind in ("truncated", "arm")
                        else seg.label if seg.kind == "marker"
                        else fs.italic(seg.label))
                ax.text((seg.start + seg.end) / 2,
                        y + sign * (height / 2 + 0.10), text,
                        ha="center", va=va, fontsize=fs.SIZE_SMALL,
                        color=fs.INK if seg.label in CASSETTE_ORDER
                        or seg.kind == "marker" else fs.INK_MUTED)

    draw_tier("construct", layout["construct"], above=True)
    draw_tier("wt", layout["wt"])
    draw_tier("mutant", layout["mutant"], gap=layout.get("break"))

    if layout.get("break"):
        _break_mark(ax, *layout["break"], _TIER_Y["mutant"], height)

    # The homology itself, bracketed on the target: it opens 15 bp ahead of the
    # start codon and closes inside the CDS, which is where the primers put it,
    # so the crossover legs land on that span and not on the whole gene.
    y_arm = _TIER_Y["wt"] + height / 2 + 0.11
    ax.plot([arm_lo, arm_hi], [y_arm, y_arm], color=fs.INK, lw=1.0,
            solid_capstyle="butt", zorder=4)
    for x in (arm_lo, arm_hi):
        ax.plot([x, x], [y_arm, y_arm - 0.07], color=fs.INK, lw=1.0, zorder=4)
    y0 = _TIER_Y["construct"] - height / 2
    for a, b in ((arm_lo, arm_hi), (arm_hi, arm_lo)):
        ax.plot([a, b], [y0, y_arm], color=fs.INK, lw=0.9, zorder=4,
                solid_capstyle="round")

    for name, text in (
            ("construct", f"{vector_name}::{_prime('rpoS')}"),
            ("wt", "HS-3 chromosome"),
            ("mutant", f"HS-3 {fs.RPOS_MINUS}")):
        # The construct tier labels its genes above the backbone, so its
        # title has to clear them rather than sit at the same height.
        lift = 0.36 if name == "construct" else 0.12
        ax.text(lo - pad, _TIER_Y[name] + height / 2 + lift, text,
                ha="left", va="bottom", fontsize=fs.SIZE_TICK, color=fs.INK)

    y = _TIER_Y["mutant"]
    complete = next(s for s in layout["mutant"]
                    if s.kind == "gene" and s.label == "rpoS")
    ax.annotate(
        "Lacks start codon",
        xy=((complete.start + complete.end) / 2, y + height / 2),
        xytext=((complete.start + complete.end) / 2, y + height / 2 + 0.38),
        ha="center", va="bottom", fontsize=fs.SIZE_SMALL, color=fs.INK,
        arrowprops=dict(arrowstyle="-", color=fs.INK_SECONDARY, lw=0.6))
    mutant_end = max(s.end for s in layout["mutant"])
    _scale_bar(ax, mutant_end - 1000, -0.66, 1000, "1 kb")


# ── The vector ────────────────────────────────────────────────
# `figure.plasmid` renders the construct as a standalone figure. The
# manuscript's panel D includes the same construct in its integration map.

def draw_plasmid(ax, size: int, features: list, *, name: str = "pRE118",
                 palette: dict = None, insert: tuple = None,
                 insert_label: str = None) -> None:
    """The sequenced vector as a ring, with the cloning site marked.

    `insert` is the (start, end) of the stretch the homology arm replaced, in
    the vector's own coordinates.
    """
    _bare(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-1.52, 1.52)
    ax.set_ylim(-1.22, 1.22)
    radius, width = 0.90, 0.085
    palette = palette or vector_palette(f.label for f in features)

    def angle(bp):
        return math.pi / 2 - 2 * math.pi * bp / size

    ax.add_patch(patches.Wedge((0, 0), radius, 0, 360, width=width,
                               facecolor=GRID_RING, edgecolor=fs.AXIS,
                               linewidth=0.5, zorder=1))
    ax.text(0, 0.07, name, ha="center", va="center", fontsize=fs.SIZE_BODY,
            fontweight="bold", color=fs.INK)
    ax.text(0, -0.09, f"{size:,} bp", ha="center", va="center",
            fontsize=fs.SIZE_TICK, color=fs.INK_SECONDARY)

    for feature in features:
        colour = palette.get(feature.label, VECTOR_FILL)
        text = (_prime(feature.label[:-1]) if feature.label.endswith("'")
                else feature.label)
        for lo, hi in feature.parts:
            ax.add_patch(patches.Wedge(
                (0, 0), radius, math.degrees(angle(hi)),
                math.degrees(angle(lo)), width=width, facecolor=colour,
                edgecolor="none", zorder=3))
        span = feature.parts[0] if len(feature.parts) > 1 else \
            (feature.parts[0][0], feature.parts[-1][1])
        _radial_label(ax, angle(sum(span) / 2), radius + 0.10, text,
                      fontsize=fs.SIZE_SMALL, color=fs.INK_SECONDARY)

    if insert:
        a = angle(sum(insert) / 2)
        r0, r1 = radius - width - 0.02, radius - width - 0.20
        ax.plot([r0 * math.cos(a), r1 * math.cos(a)],
                [r0 * math.sin(a), r1 * math.sin(a)],
                color=fs.INK, lw=0.9, zorder=4)
        ax.add_patch(patches.Wedge((0, 0), radius, math.degrees(angle(insert[1])),
                                   math.degrees(angle(insert[0])), width=width,
                                   facecolor=fs.INK, edgecolor="none", zorder=4))
        _radial_label(ax, a, r1 - 0.02, insert_label or "cloning site",
                      inward=True, fontsize=fs.SIZE_SMALL, color=fs.INK)
