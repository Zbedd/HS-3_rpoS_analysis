"""Panels for the rpoS-locus figure. Each draws into an axes the caller owns.

No function here creates a figure, writes a file, or sets an rcParam —
`viz` owns all of that, and the caller owns the page.

Linear gene tracks are drawn in base pairs, with omitted sequence marked by
an axis break. Circular maps show feature positions around the replicon.
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


def _vector_label(label: str) -> str:
    """Format gene symbols while preserving annotation keys for colours."""
    if label in ("KanR", "kanR", "sacB"):
        return fs.italic("kanR" if label == "KanR" else label)
    return label


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


def _defective_rpos(ax, seg, y, height, *, terminus):
    """A broken terminal edge marks a partial copy without changing its span."""
    span = seg.end - seg.start
    notch = min(80, span * 0.16)
    half = height / 2
    if terminus == "C":
        points = [(0, -half), (span, -half), (span - notch, -half / 2),
                  (span, 0), (span - notch, half / 2), (span, half), (0, half)]
    else:
        head = min(240, span * 0.35)
        points = [(0, -half), (span - head, -half), (span, 0),
                  (span - head, half), (0, half), (notch, half / 2),
                  (0, 0), (notch, -half / 2)]
    vertices = [(seg.start + x if seg.strand > 0 else seg.end - x, y + dy)
                for x, dy in points]
    ax.add_patch(patches.Polygon(vertices, closed=True,
                                 facecolor=cassette_fill("rpoS"),
                                 edgecolor=fs.INK, hatch="///",
                                 linewidth=0.5, zorder=2))


def draw_integration(ax, layout: dict, *, palette: dict = None,
                     vector_name: str = "pRE118") -> None:
    """Design over outcome: the construct, the target, and what resulted.

    The circular construct sits above two chromosome tracks. Those tracks
    share a base-pair scale outside the marked break in the integrant.
    """
    _bare(ax)
    height = 0.34
    arm_lo, arm_hi = layout["arm"]
    palette = palette or {}

    lo = min(s.start for tier in ("wt", "mutant")
             for s in layout[tier])
    hi = max(s.end for tier in ("wt", "mutant") for s in layout[tier])
    pad = 0.03 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-0.92, 4.20)

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
                _arrow(ax, seg.start, seg.end, seg.strand, y, height,
                       palette.get(seg.label, VECTOR_FILL), edgecolor=fs.INK,
                       min_head=240)
            elif seg.kind == "truncated":
                _defective_rpos(ax, seg, y, height, terminus="C")
            elif name == "mutant" and seg.kind == "gene" and seg.label == "rpoS":
                _defective_rpos(ax, seg, y, height, terminus="N")
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
                        else _vector_label(seg.label) if seg.kind == "marker"
                        else fs.italic(seg.label))
                ax.text((seg.start + seg.end) / 2,
                        y + sign * (height / 2 + 0.10), text,
                        ha="center", va=va, fontsize=fs.SIZE_SMALL,
                        color=fs.INK if seg.label in CASSETTE_ORDER
                        or seg.kind == "marker" else fs.INK_MUTED)

    size, features = layout["circular_construct"]
    shown = [f for f in features
             if f.label in ("rpoS'", "KanR", "sacB")]
    arm_feature = next(f for f in shown if f.label == "rpoS'")
    rotation = -math.pi / 2 - 2 * math.pi * (arm_feature.start + arm_feature.end) / (2 * size)
    centre = ((arm_lo + arm_hi) / 2 - (lo - pad)) / (hi - lo + 2 * pad)
    plasmid_ax = ax.inset_axes([centre - 0.19, 0.56, 0.38, 0.43])
    draw_plasmid(plasmid_ax, size, shown, name=vector_name,
                 palette={**palette, "rpoS'": cassette_fill("rpoS")},
                 rotation=rotation, clockwise=False,
                 arm_inside=True, directional=True)
    draw_tier("wt", layout["wt"])
    draw_tier("mutant", layout["mutant"], gap=layout.get("break"))

    if layout.get("break"):
        _break_mark(ax, *layout["break"], _TIER_Y["mutant"], height)

    # Follow the crossed loop from the left chromosome into the vector after
    # the arm, through kanR then sacB, and back to the right chromosome.
    y_arm = _TIER_Y["wt"] + height / 2 + 0.11
    ax.plot([arm_lo, arm_hi], [y_arm, y_arm], color=fs.INK, lw=1.0,
            solid_capstyle="butt", zorder=4)
    for x in (arm_lo, arm_hi):
        ax.plot([x, x], [y_arm, y_arm - 0.07], color=fs.INK, lw=1.0, zorder=4)
    for bp, target in ((arm_feature.end, arm_lo), (arm_feature.start, arm_hi)):
        angle = rotation + 2 * math.pi * bp / size
        ax.add_artist(patches.ConnectionPatch(
            xyA=(0.90 * math.cos(angle), 0.90 * math.sin(angle)),
            xyB=(target, y_arm), coordsA="data", coordsB="data",
            axesA=plasmid_ax, axesB=ax, color=fs.INK, lw=0.9, zorder=4))

    for name, text in (
            ("wt", "HS-3 chromosome"),
            ("mutant", f"HS-3 {fs.RPOS_MINUS}")):
        ax.text(lo - pad, _TIER_Y[name] + height / 2 + 0.12, text,
                ha="left", va="bottom", fontsize=fs.SIZE_TICK, color=fs.INK)

    y = _TIER_Y["mutant"]
    for kind, text, terminus in (("truncated", "C-terminal truncation", "C"),
                                 ("gene", "Lacks start codon", "N")):
        copy = next(s for s in layout["mutant"]
                    if s.kind == kind and s.label == "rpoS")
        edge = copy.end if (terminus == "C") == (copy.strand > 0) else copy.start
        ax.annotate(
            text, xy=(edge, y + height / 2),
            xytext=((copy.start + copy.end) / 2, y + height / 2 + 0.24),
            ha="center", va="bottom", fontsize=fs.SIZE_SMALL, color=fs.INK,
            arrowprops=dict(arrowstyle="-", color=fs.INK_SECONDARY, lw=0.6))
    mutant_end = max(s.end for s in layout["mutant"])
    _scale_bar(ax, mutant_end - 1000, -0.66, 1000, "1 kb")


# ── The vector ────────────────────────────────────────────────
# `figure.plasmid` renders the construct as a standalone figure. The
# manuscript's panel D includes the same construct in its integration map.

def draw_plasmid(ax, size: int, features: list, *, name: str = "pRE118",
                 palette: dict = None, insert: tuple = None,
                 insert_label: str = None, rotation: float = math.pi / 2,
                 arm_inside: bool = False, directional: bool = False,
                 clockwise: bool = True) -> None:
    """The sequenced vector as a ring, with the cloning site marked.

    `insert` is the (start, end) of the stretch the homology arm replaced, in
    the vector's own coordinates.
    """
    _bare(ax)
    ax.set_aspect("equal")
    ax.set_xlim(-1.52, 1.52)
    ax.set_ylim(-1.22, 1.22)
    radius, width = 0.90, 0.13 if directional else 0.085
    palette = palette or vector_palette(f.label for f in features)

    def angle(bp):
        return rotation + (-1 if clockwise else 1) * 2 * math.pi * bp / size

    def arc(lo, hi):
        return sorted((math.degrees(angle(lo)), math.degrees(angle(hi))))

    ax.add_patch(patches.Wedge((0, 0), radius, 0, 360, width=width,
                               facecolor=GRID_RING, edgecolor=fs.AXIS,
                               linewidth=0.5, zorder=1))
    ax.text(0, 0.18 if arm_inside else 0.07, name,
            ha="center", va="center", fontsize=fs.SIZE_BODY,
            fontweight="bold", color=fs.INK)
    ax.text(0, -0.16 if arm_inside else -0.09, f"{size:,} bp",
            ha="center", va="center",
            fontsize=fs.SIZE_TICK, color=fs.INK_SECONDARY)

    for feature in features:
        colour = palette.get(feature.label, VECTOR_FILL)
        text = (_prime(feature.label[:-1]) if feature.label.endswith("'")
                else _vector_label(feature.label))
        for lo, hi in feature.parts:
            if directional and feature.label in ("KanR", "sacB", "traJ"):
                start, end = angle(lo), angle(hi)
                tail, tip = (start, end) if feature.strand > 0 else (end, start)
                direction = 1 if tip > tail else -1
                head = min(abs(tip - tail) * 0.35, 0.20)
                shoulder = tip - direction * head
                angles = [tail + (shoulder - tail) * i / 40 for i in range(41)]
                vertices = [(radius * math.cos(a), radius * math.sin(a)) for a in angles]
                vertices.append(((radius - width / 2) * math.cos(tip),
                                 (radius - width / 2) * math.sin(tip)))
                vertices.extend(((radius - width) * math.cos(a),
                                 (radius - width) * math.sin(a)) for a in reversed(angles))
                ax.add_patch(patches.Polygon(vertices, closed=True,
                                             facecolor=colour, edgecolor="none", zorder=3))
                continue
            ax.add_patch(patches.Wedge(
                (0, 0), radius, *arc(lo, hi), width=width, facecolor=colour,
                edgecolor="none", zorder=3))
        span = feature.parts[0] if len(feature.parts) > 1 else \
            (feature.parts[0][0], feature.parts[-1][1])
        if arm_inside and feature.label == "rpoS'":
            ax.text(0, -0.58, text, ha="center", va="center",
                    fontsize=fs.SIZE_SMALL, color=fs.INK)
            continue
        _radial_label(ax, angle(sum(span) / 2), radius + 0.10, text,
                      fontsize=fs.SIZE_SMALL, color=fs.INK_SECONDARY)

    if insert:
        a = angle(sum(insert) / 2)
        r0, r1 = radius - width - 0.02, radius - width - 0.20
        ax.plot([r0 * math.cos(a), r1 * math.cos(a)],
                [r0 * math.sin(a), r1 * math.sin(a)],
                color=fs.INK, lw=0.9, zorder=4)
        ax.add_patch(patches.Wedge((0, 0), radius, *arc(*insert), width=width,
                                   facecolor=fs.INK, edgecolor="none", zorder=4))
        _radial_label(ax, a, r1 - 0.02, insert_label or "cloning site",
                      inward=True, fontsize=fs.SIZE_SMALL, color=fs.INK)
