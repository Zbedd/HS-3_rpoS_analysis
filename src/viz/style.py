"""HS-3 house figure style — one style, reused by every figure in this repo.

Every plotting module (CR, SURV) imports this module and takes
its typography, colours, marks, geometry and export settings from here. Change
a value here and every figure changes with it.

Panel arrangement lives next door in `viz.layout`; both are re-exported
by `viz`, which is what plotting modules import.

    import viz as fs

    fs.use()                                  # install the rcParams, once
    fig, ax = plt.subplots(figsize=(fs.COL_DOUBLE, 3.0))
    ax.plot(x, y, color=fs.colour("rpoS"), marker=fs.marker("rpoS"),
            label=fs.label("rpoS"))
    fs.finish_axis(ax)
    fs.legend_outside(ax, title="Strain")     # legends sit beside the data
    fs.save(fig, path)                        # PNG (600 dpi) + vector PDF


CONVENTIONS
-----------
Geometry.  Figures are sized in journal column widths (COL_SINGLE /
    COL_ONEHALF / COL_DOUBLE), so type lands at its real printed size.

Type.  One family (Arial, with Helvetica/DejaVu fallbacks), four sizes
    (SIZE_TITLE / SIZE_BODY / SIZE_TICK / SIZE_SMALL) derived from BASE_PT.
    Nothing sets a size off that ladder. Titles and axis labels share
    SIZE_BODY and separate by weight rather than by half a point. Gene
    symbols go through `italic()` / `label()` rather than hand-written
    mathtext. PDF/PS text is embedded as TrueType (fonttype 42), so vector
    output stays editable in Illustrator/Inkscape.

Colour.  Series colour is identity, keyed on the strain, from `COLOURS`. The
    same strain is the same colour in every assay; colour does not encode
    magnitude and is not assigned by rank or list position. Values come from
    PALETTE PROVENANCE below.

Marks.  Every strain also carries a distinct marker shape (`MARKERS`), so
    identity survives greyscale and colour-vision deficiency. Three palette
    slots sit below 3:1 contrast on white, which is legal only when identity
    is also carried by shape or a direct axis label. The marker is the datum
    and the line only joins them, so the marker is the heavier mark of the
    two (`lines.markersize` against `lines.linewidth` / `ERRORBAR_LW`).

Chrome.  White surface, no top/right spines, a hairline horizontal grid behind
    the data, chrome greys lighter than any text.

Text.  Two inks, and only two: INK names things (titles, axis labels, panel
    letters, legend titles), INK_SECONDARY reads them out (tick labels,
    legend entries, annotations, footnotes). INK_MUTED is a mark colour, not
    a text colour.

Legends.  Never over the data. `legend()` places the key in a corner the data
    leaves empty, which is where a panel with a free corner should carry it;
    `legend_outside()` gives it its own column to the right when no corner is
    free, and both layout engines reserve that column inside the figure rather
    than growing its width.

Export.  `save()` writes a 600 dpi PNG and a vector PDF from one call, with
    identical margins.


PALETTE PROVENANCE
------------------
The categorical hues carry the hue angles of the house swatch set in
`PASTELS`, re-rendered at series strength. A swatch cannot serve as a series
colour at its own value: the set runs chroma 0.05-0.10 against the 0.10 floor a
categorical palette needs, three swatches sit above the 0.43-0.77 lightness
band, and two of the ten (pale yellow / target beige) are the same colour to
within dE 2. So each slot keeps its swatch's hue and takes the lightness and
chroma nearest that swatch which still clears both hard gates:

    CVD separation      >= 9.2   (OKLab dE x100, min of protanopia/deuteranopia
                                  under Machado-Oliveira-Fernandes 2009 @ 1.0)
    normal-vision floor >= 15.6  (same metric, unsimulated)

Those two thresholds are the margins the previous palette reached, so the swap
cannot cost separation. They are applied over **all pairs** of each sequence,
not just adjacent ones: a reader compares any two lines in a survival or CFU
plot. The previous palette passed only adjacent pairs; a non-adjacent pair
had a CVD separation of 3.2.

Two slots are pinned to the swatch exactly, because the CR figures are
published on them: `orange` #FE8C26 and `blue` #4A5FAA. `blue` therefore
belongs to PA14 rather than WT, so that CR's lineage pair is those two hexes;
WT takes `violet`. Verified with the skill's validator on the sequence each
figure really draws, all pairs:

    SURV lines    WT,rpoS- (genotype keying)     CVD 30.2 / normal 37.3
    CR (lineage)  HS-3,PA14                      CVD 30.2 / normal 37.3

CR draws by lineage, not by strain — `panels` colours every mark through
`lineage_colour` and hands genotype to `fill_for` — so CR needs one separated
pair, not five separated hues.

Three slots WARN on contrast against white (aqua 2.25:1, orange 2.34:1,
violet 2.89:1); `MARKERS` and the direct group labels on every categorical
axis are the secondary encoding that covers it. No slot is weaker than the
2.17:1 the previous palette accepted.

The swatch set has seven distinct hue angles and no green, so `green` is the
one slot with no swatch behind it; it fills the 100-160 deg gap. That also
leaves no hue free for `ORDINAL_FAMILIES`, whose teal ramp now shares a hue
family with `aqua` — see the note there.

To add a strain, give it a hue angle from `PASTELS`, re-render it into the
lightness band, and re-run the validator on the sequences the figures draw.

    python <dataviz-skill>/scripts/validate_palette.py \
        "#9291d8,#FE8C26,#4dbeb3,#186c24" --mode light --surface "#ffffff" \
        --pairs all
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402


# The public surface, re-exported by `viz`. A name absent here is
# internal to this module.
__all__ = [
    "COL_SINGLE", "COL_ONEHALF", "COL_DOUBLE", "ROW_GAP", "PANEL_GAP",
    "LETTER_BAND", "INK", "INK_SECONDARY",
    "INK_MUTED", "AXIS", "GRID", "SURFACE", "PASTELS", "HUES", "SERIES_CYCLE",
    "COLOURS", "ALIASES", "MARKERS", "LINESTYLES", "LINEAGE",
    "LINEAGE_COLOURS", "lineage", "lineage_colour", "MUTANTS",
    "is_mutant", "fill_for", "COLORMAPS", "colormap", "DOMAIN_COLOURS", "DOMAIN_UNMAPPED", "domain_colour",
    "ROLES", "role", "role_edge", "ORDINAL_FAMILIES", "ORDINAL_ENDS",
    "ordinal", "ramp_ends", "SIGMA70_ORDER",
    "GENOTYPE_COLOURS", "genotype_colour", "RPOS_MINUS",
    "ink_on", "italic", "lesion", "LABELS", "colour", "marker", "linestyle",
    "label",
    "arrow", "scale_bar",
    "BASE_PT", "SIZE_TITLE", "SIZE_BODY", "SIZE_TICK", "SIZE_SMALL",
    "ERRORBAR_LW", "RC_PARAMS", "use", "finish_axis", "fit_xticklabels",
    "legend", "legend_outside", "footnote", "reference_line",
    "SAVE_FORMATS", "save",
]


# ── Geometry ──────────────────────────────────────────────────
# Journal column widths in inches. Draw at the width the figure will be
# printed at; never draw huge and let the typesetter shrink the type.
COL_SINGLE = 3.46    #  88 mm
COL_ONEHALF = 4.72   # 120 mm
COL_DOUBLE = 7.09    # 180 mm

# What separates the panels of a stacked figure, in inches. A figure states
# its panels' heights and takes these; nothing hand-tunes a gap.
ROW_GAP = 0.06       # between two rows of one panel
PANEL_GAP = 0.26     # clear space under one panel before the next begins
LETTER_BAND = 0.24   # room above a panel for its letter


# ── Ink & chrome ──────────────────────────────────────────────
# Two text inks, two chrome greys, one neutral mark colour. Every string in
# every figure is INK or INK_SECONDARY; nothing else is a text colour.
INK = "#0b0b0b"            # titles, axis labels, panel letters, legend titles
INK_SECONDARY = "#52514e"  # tick labels, legend entries, annotations, footnotes
INK_MUTED = "#898781"      # marks only: the neutral series colour, reference rules
AXIS = "#c3c2b7"           # spines and tick marks
GRID = "#e1e0d9"           # hairline grid behind the data
SURFACE = "#ffffff"        # figure + axes background


# ── House swatches ────────────────────────────────────────────
# The set this repository's palette is built from. These are *fill* values —
# large labelled areas, cartoon blocks, a wash behind a span. They are too pale
# and too weak in chroma to mark a series (see PALETTE PROVENANCE), so a line,
# a marker or a violin takes its colour from `HUES` instead, which carries
# these hue angles at series strength.
#
# `beige` and `pale_yellow` are the same hue to within dE 2; `periwinkle`,
# `medium_blue` and `dark_blue` are one hue at three lightnesses, which is why
# the set supplies seven hue angles, not ten.
PASTELS = {
    "pale_yellow": "#F5DDB6",
    "peach":       "#F5CBBC",
    "pink":        "#F5B6BF",
    "mauve":       "#CEA1B6",
    "periwinkle":  "#9B9CC0",
    "medium_blue": "#728AB9",
    "dark_blue":   "#4A5FAA",
    "teal":        "#4FB8AE",
    "beige":       "#FEDFAB",
    "orange":      "#FE8C26",
}


# ── Categorical palette ───────────────────────────────────────
# The eight series hues: each `PASTELS` hue angle rendered at the lightness and
# chroma nearest its swatch that still clears the gates in PALETTE PROVENANCE.
# `orange` and `blue` are their swatches unchanged — the CR figures are
# published on those two values. `COLOURS` picks from these; unmapped series
# fall back to this cycle.
HUES = {
    "blue":    "#4A5FAA",   # dark_blue, exact
    "orange":  "#FE8C26",   # orange, exact
    "aqua":    "#4dbeb3",   # teal
    "yellow":  "#ab7a00",   # pale_yellow / beige
    "magenta": "#9e3971",   # mauve
    "green":   "#186c24",   # no swatch: the set's one hue gap
    "violet":  "#9291d8",   # periwinkle
    "red":     "#e25878",   # pink
}
SERIES_CYCLE = list(HUES.values())

# Strain -> colour. One strain, one colour, everywhere in the repository.
# Keys are the raw sample names used in data/ and in the assay configs; the
# aliases below map each assay's local spelling onto the same identity.
COLOURS = {
    # PA14 holds `blue` and WT holds `violet`, so that CR's lineage pair is
    # the two hexes its figures are published on.
    "WT":         HUES["violet"],
    "kanR":       HUES["orange"],
    "rpoS":       HUES["aqua"],
    "PA14":       HUES["blue"],
    "PA14_rpoS":  HUES["magenta"],
    # Media-only / blank references are achromatic: they are a baseline, not
    # a competing identity.
    "ctrl":       INK_MUTED,
    "blank":      INK_MUTED,
}

# Assay-local spellings that denote the same strain as a key above.
ALIASES = {
    "rpoS::ΩkanR": "rpoS",
    "ΔrpoS": "rpoS",
    "PA14 ΔrpoS": "PA14_rpoS",
    "PA14_ΔrpoS": "PA14_rpoS",
    "control": "ctrl",
}

# Marker shape per strain — the secondary encoding for greyscale, CVD, and the
# sub-3:1 palette slots.
MARKERS = {
    "WT":         "o",
    "kanR":       "s",
    "rpoS":       "^",
    "PA14":       "v",
    "PA14_rpoS":  "P",
    "ctrl":       "X",
    "blank":      "X",
}

LINESTYLES = {}

# ── Lineage palette ───────────────────────────────────────────
# One colour per organism, for a figure that pairs a control with its mutant
# and lets fill carry the genotype. Each organism takes the hue its reference
# strain already has in COLOURS, so the control keeps its identity and only
# the mutant is recoloured.
#
#   orange/blue, to each other   CVD 30.2 · normal 37.3 — both gates pass
#
# A figure keys on this or on COLOURS, never both: per-strain colour separates
# four strains, lineage colour separates two organisms and hands genotype to
# `fill_for`.
LINEAGE = {
    "WT": "HS-3", "kanR": "HS-3", "rpoS": "HS-3",
    "PA14": "PA14", "PA14_rpoS": "PA14",
}
LINEAGE_COLOURS = {"HS-3": HUES["orange"], "PA14": HUES["blue"]}


def lineage(name: str) -> str:
    """Organism a strain belongs to."""
    return LINEAGE.get(_canonical(name), _canonical(name))


def lineage_colour(name: str, default: str = INK_SECONDARY) -> str:
    """Colour of the organism a strain belongs to."""
    return LINEAGE_COLOURS.get(lineage(name), default)


# ── Genotype fill ─────────────────────────────────────────────
# A second channel over the strain colour, for a figure that pairs a control
# with its mutant: both keep the strain's own hue, and only the face changes.
# Colour still carries identity, so a figure may use this or ignore it.
MUTANTS = frozenset({"rpoS", "PA14_rpoS"})


def is_mutant(name: str) -> bool:
    """Whether a strain carries a lesion."""
    return _canonical(name) in MUTANTS


# A third keying, for a figure drawing one organism's two genotypes that cannot
# spend fill on the genotype — SURV's open marker already means a timepoint
# containing a non-detect. Control and mutant take the two hues the CR figures
# are published on, so the two assays' main panels read as a pair.
#
#   orange/blue, to each other   CVD 30.2 · normal 37.3 — both gates pass
#
# A figure keys on this, on COLOURS, or on LINEAGE_COLOURS — never two at once.
GENOTYPE_COLOURS = {"control": HUES["orange"], "mutant": HUES["blue"]}


def genotype_colour(name: str) -> str:
    """Colour separating a control from its mutant, ignoring which strain."""
    return GENOTYPE_COLOURS["mutant" if is_mutant(name) else "control"]


def fill_for(name: str, tone: str = None) -> tuple:
    """(face, edge) drawing `name` open when unmutated, filled when mutant.

    The edge is the colour either way, so the two states read as one series.
    `tone` overrides the colour, which is otherwise the strain's own.
    """
    c = tone if tone is not None else colour(name)
    return (c if is_mutant(name) else SURFACE), c


# ── Continuous maps ───────────────────────────────────────────
# For image fields rather than categorical marks, registered by name so a
# figure picks one instead of hardcoding a colormap.
#
# `od_chrom` diverges about zero because zero is the agar reference the metric
# is built on, and 28% of crop pixels sit below it. RdBu_r is ColorBrewer's
# colour-vision-safe diverging scheme: it gives up the perceptual uniformity of
# a sequential map to keep a readable zero, which a referenced field needs.
# Positive — colony redder than its own agar — runs red, matching the assay.
COLORMAPS = {"od_chrom": "RdBu_r"}


def colormap(name: str):
    """Colormap registered under `name`."""
    return plt.get_cmap(COLORMAPS[name])


# ── Protein-domain palette ────────────────────────────────────
# For the domain-architecture cartoons. Keyed on Pfam accession, so the same
# domain is the same colour in every figure that draws it — RQ1 panel B and
#
# Validated on a harder test than the series palette: a cartoon shows whichever
# domains a protein happens to carry, so the check is ALL pairs, not adjacent
# ones. Within each family the hues follow the domain's N→C position, so colour
# also tracks sequence order.
#
#   σ⁷⁰      (6 domains, all pairs)   CVD  6.9 · normal 15.6  — normal-vision
#       gate passes; CVD sits in the 6–8 floor band, which is legal only with
#       secondary encoding. The in-block domain labels are that encoding.
#
# Orange is absent from both sets: it collapses against magenta and green under
# deuteranopia, and no six-hue σ⁷⁰ set containing it clears the normal-vision
# floor. Re-run the search in the module header before adding a domain.
DOMAIN_COLOURS = {
    # σ⁷⁰ regions, N→C
    "PF03979": HUES["blue"],     # Sigma70_r1_1
    "PF00140": HUES["aqua"],     # Sigma70_r1_2
    "PF04542": HUES["yellow"],   # Sigma70_r2
    "PF04539": HUES["green"],    # Sigma70_r3
    "PF04545": HUES["violet"],   # Sigma70_r4
    "PF08281": HUES["red"],      # Sigma70_r4_2
}

# Any domain outside the two families above: neutral, so an unclassified block
# does not read as an identity.
DOMAIN_UNMAPPED = "#cfcec7"


# σ⁷⁰ regions in N→C order. They are a sequence, not a set: the order is the
# information, and a lightness ramp carries it where categorical hues throw it
# away. Taking them off the categorical palette also returns four hues to
# strain duty — see the roles block below.
SIGMA70_ORDER = ("PF03979", "PF00140", "PF04542", "PF04539", "PF04545",
                 "PF08281")

# The non-essential region: an insertion between r1.2 and r2 that Group 1
# factors carry and Group 2 factors lack, so it appears on RpoD and on nothing
# else the figures draw. It takes the ramp value midway between the regions it
# separates, which joins it to the series without moving their six steps.
SIGMA70_NCR = "PF04546"

# σ54 domains, N->C. RpoN shares no fold with the σ70 family — a different
# protein with three domains of its own — so it is its own ordered series
# rather than extra steps on the σ70 ramp. The two never appear in one panel.
SIGMA54_ORDER = ("PF00309", "PF04963", "PF04552")


def domain_colour(accession: str) -> str:
    """Block colour for a Pfam accession in a domain-architecture cartoon.

    Ordered families resolve to a ramp step; everything else keeps its
    categorical hue.
    """
    for order in (SIGMA70_ORDER, SIGMA54_ORDER):
        if accession in order:
            return ordinal(order.index(accession), len(order))
    if accession == SIGMA70_NCR:
        anchor = SIGMA70_ORDER.index("PF00140")
        return ordinal(2 * anchor + 1, 2 * len(SIGMA70_ORDER) - 1)
    return DOMAIN_COLOURS.get(accession, DOMAIN_UNMAPPED)


# ── Roles ─────────────────────────────────────────────────────
# Within one figure, a hue means exactly one thing. That is the rule colour
# obeys here — not the absence of hue. What went wrong before was narrower and
# fixable: rpoS was aqua on one panel and violet on the next, adjacent panels
# reused a hue for unrelated things, and vector blocks took whatever hue their
# index landed on.
#
# So hue still carries identity. These neutrals are for what identity is *not*
# being claimed about — context a panel draws for orientation, and the span of
# a replicon nothing is annotated on.
#
#   backdrop  drawn, but not part of the claim — context genes
#   feature   drawn and named, but not the subject
#   focus     a subject drawn deliberately without a hue
#   ring      the span of a circular replicon nothing is annotated on
ROLES = {
    "backdrop": "#eae8e0",
    "feature":  "#c0bcae",
    "focus":    "#6b6759",
    "ring":     "#f6f5f1",
}
ROLE_EDGES = {
    "backdrop": INK_MUTED,
    "feature":  INK,
    "focus":    INK,
    "ring":     AXIS,
}


def role(name: str) -> str:
    """Fill for a drawn feature, by the part it plays rather than what it is."""
    return ROLES[name]


def role_edge(name: str) -> str:
    """Outline for a drawn feature in the given role."""
    return ROLE_EDGES[name]


# ── Colour space ──────────────────────────────────────────────
# OKLab, for deriving one colour from another instead of picking a second one
# by eye. Only `ramp_ends` needs this; everything else in the module names its
# colours outright.

def _srgb_to_linear(c: float) -> float:
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _linear_to_srgb(c: float) -> float:
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def _to_oklab(hex_colour: str) -> tuple:
    h = hex_colour.lstrip("#")
    r, g, b = (_srgb_to_linear(int(h[i:i + 2], 16) / 255) for i in (0, 2, 4))
    l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
    m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
    s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
    l_, m_, s_ = (v ** (1 / 3) if v >= 0 else -((-v) ** (1 / 3))
                  for v in (l, m, s))
    return (0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
            1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
            0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_)


def _from_oklab(lightness: float, a: float, b: float) -> tuple:
    l_ = lightness + 0.3963377774 * a + 0.2158037573 * b
    m_ = lightness - 0.1055613458 * a - 0.0638541728 * b
    s_ = lightness - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    return (4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
            -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
            -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s)


def _polar(hex_colour: str) -> tuple:
    """(lightness, chroma, hue) of a colour, hue in radians."""
    lightness, a, b = _to_oklab(hex_colour)
    return lightness, math.hypot(a, b), math.atan2(b, a)


def _hex_at(lightness: float, chroma: float, hue: float) -> tuple:
    """(hex, whether it left the sRGB gamut) at this lightness/chroma/hue."""
    rgb = _from_oklab(lightness, chroma * math.cos(hue), chroma * math.sin(hue))
    clipped = any(v < -1e-4 or v > 1 + 1e-4 for v in rgb)
    return "#%02x%02x%02x" % tuple(
        max(0, min(255, round(_linear_to_srgb(max(0.0, min(1.0, v))) * 255)))
        for v in rgb), clipped


def ramp_ends(anchor: str, dark_lightness: float = 0.28,
              chroma_gain: float = 1.15) -> tuple:
    """(anchor, a dark end in the anchor's own hue) — a ramp for `ordinal`.

    The light end *is* the anchor, so an ordered series drawn on this ramp
    reads as belonging to whatever the anchor identifies. The dark end is
    derived, not chosen: same hue, `dark_lightness`, and as much chroma as the
    gamut allows up to `chroma_gain` x the anchor's, which keeps the deep steps
    from washing out to grey.

    Defaults preserve ordered-domain contrast (monotone lightness,
    adjacent dL >= 0.06, light end >= 2:1 on white, one hue).
    """
    lightness, chroma, hue = _polar(anchor)
    ceiling = chroma * chroma_gain
    low, high = 0.0, ceiling
    for _ in range(32):                     # largest in-gamut chroma <= ceiling
        mid = (low + high) / 2
        if _hex_at(dark_lightness, mid, hue)[1]:
            high = mid
        else:
            low = mid
    return anchor, _hex_at(dark_lightness, low, hue)[0]


# Light to dark along an ordered series. A ramp reads as an order, which a
# categorical set cannot, and because it is read as an order it does not
# compete with the categorical hues a neighbouring panel is using.
ORDINAL_FAMILIES = {
    # Protein domains N->C, on rpoS's own colour. RpoS's domains are parts of
    # the thing `rpoS` names, so they ramp in its hue: the gene arrow on a
    # locus map and the domain blocks on a cartoon read as one family.
    # Derived from COLOURS, so recolouring the strain carries the ramp with it
    # and neither can drift from the other.
    #
    # This is also what keeps the ramp off the other seven hues, which it
    # clears by >= 9.3 (OKLab dE x100). It may not be reused for an ordered
    # series that is not RpoS — that would put rpoS's hue on something else.
    "domain": ramp_ends(COLOURS["rpoS"]),
}
ORDINAL_ENDS = ORDINAL_FAMILIES["domain"]


def ordinal(index: int, count: int, ends=None, family: str = "domain") -> str:
    """Step `index` of `count` along an ordered ramp, light to dark.

    For a variable whose order is the point — a protein's domains N to C, a
    series of timepoints. A reader gets the ordering for free, and it survives
    every kind of colour-vision deficiency because it is carried by lightness.
    """
    low, high = ends or ORDINAL_FAMILIES[family]
    if count <= 1:
        return high
    fraction = index / (count - 1)
    start = [int(low.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4)]
    stop = [int(high.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4)]
    return "#%02x%02x%02x" % tuple(
        round(a + (b - a) * fraction) for a, b in zip(start, stop))


# How much better dark ink must be before a label gives up white. Several
# palette slots are near-ties on raw contrast (blue is 1.01), so white is kept
# unless dark ink is clearly ahead.
_DARK_INK_MARGIN = 1.5


def ink_on(background: str) -> str:
    """Text colour that stays legible on `background`.

    White unless dark ink beats it by `_DARK_INK_MARGIN`. Every slot that
    keeps white clears 3.9:1 against it, above the 3:1 large-text bar for the
    bold in-block labels this serves.
    """
    red, green, blue = (int(background.lstrip("#")[i:i + 2], 16) / 255
                        for i in (0, 2, 4))

    def channel(component):
        return (component / 12.92 if component <= 0.04045
                else ((component + 0.055) / 1.055) ** 2.4)

    luminance = (0.2126 * channel(red) + 0.7152 * channel(green)
                 + 0.0722 * channel(blue))
    on_white = 1.05 / (luminance + 0.05)
    on_dark = (luminance + 0.05) / (0.0075 + 0.05)   # against INK, not pure black
    return INK if on_dark > _DARK_INK_MARGIN * on_white else "#ffffff"


# ── Marks ─────────────────────────────────────────────────────

def arrow(ax, start, end, strand, y, height, facecolor, *, edgecolor=None,
          lw=0.5, head_frac=0.35, min_head=0.0, zorder=2, **kwargs):
    """A pentagon arrow from `start` to `end`, pointing along `strand`.

    The shape every feature drawn on a sequence uses — a gene on a chromosome,
    a domain on a protein — so the same silhouette means the same thing in
    every panel. The head is a fraction of the feature's own length, so a
    short one keeps a readable body, and `min_head` caps it in data units so a
    long one does not grow a spear.
    """
    import matplotlib.patches as patches

    span = end - start
    head = min(span * head_frac, min_head) if min_head else span * head_frac
    if strand >= 0:
        body, tip = start, end
        shoulder = tip - head
    else:
        body, tip = end, start
        shoulder = tip + head
    ax.add_patch(patches.Polygon(
        [(body, y - height / 2), (shoulder, y - height / 2), (tip, y),
         (shoulder, y + height / 2), (body, y + height / 2)],
        closed=True, facecolor=facecolor, edgecolor=edgecolor or INK,
        linewidth=lw, zorder=zorder, joinstyle="miter", **kwargs))


def scale_bar(ax, x0, y, span, label, *, tick=0.06):
    """A bar `span` data units long, ticked at both ends and named beneath.

    What a cartoon panel carries instead of an axis: the reader needs the
    scale, not a coordinate, and an axis on a diagram invites reading
    positions off it that the diagram does not claim.
    """
    ax.plot([x0, x0 + span], [y, y], color=INK_SECONDARY, lw=1.0,
            solid_capstyle="butt", zorder=5)
    for x in (x0, x0 + span):
        ax.plot([x, x], [y - tick, y + tick], color=INK_SECONDARY, lw=1.0,
                zorder=5)
    ax.text(x0 + span / 2, y - tick - 0.04, label, ha="center", va="top",
            fontsize=SIZE_SMALL, color=INK_SECONDARY)


# ── Display labels ────────────────────────────────────────────

def italic(text: str) -> str:
    """Render `text` in the italic of the figure font, via mathtext.

    Use for gene symbols. `\\mathit` on a custom mathtext fontset picks up
    `font.sans-serif`, so italics match the surrounding type rather than
    dropping into Computer Modern.
    """
    return rf"$\mathit{{{text}}}$"


def lesion(gene: str) -> str:
    """Gene symbol carrying a lesion: italic symbol, superscript minus."""
    return rf"$\mathit{{{gene}}}^{{-}}$"


RPOS_MINUS = lesion("rpoS")

# STRAIN NAMING — the one place a strain gets the name a reader sees.
#
# Keys are the raw sample names in data/ and in the assay configs. They are
# the experiment's names and are left alone; the values are the paper's, and
# the two deliberately disagree:
#
#   data `WT`     → "Parental"   the unmarked isolate, no cassette
#   data `kanR`   → "WT"         the kanamycin-marked strain everything is
#                                compared against, which the paper treats as
#                                wild type
#   data `rpoS`   → rpoS⁻        the rpoS mutant
#
# So "WT" in a figure is the marked strain, and the unmarked one it descends
# from is "Parental". The rename is display-only: colours, markers and every
# key in COLOURS above still use the raw names, so no figure changes hue.
LABELS = {
    "WT":         "Parental",
    "kanR":       "WT",
    "rpoS":       RPOS_MINUS,
    "PA14":       "WT",
    "PA14_rpoS":  RPOS_MINUS,
    "ctrl":       "ctrl",
    "blank":      "blank",
}


def _canonical(name: str) -> str:
    return ALIASES.get(name, name)


def colour(name: str, default: str = INK_SECONDARY) -> str:
    """Series colour for a strain (accepts any assay's local spelling)."""
    return COLOURS.get(_canonical(name), default)


def marker(name: str, default: str = "o") -> str:
    """Marker shape for a strain."""
    return MARKERS.get(_canonical(name), default)


def linestyle(name: str, default: str = "-"):
    """Line style for a strain."""
    return LINESTYLES.get(_canonical(name), default)


def label(name: str) -> str:
    """What a figure calls a strain — see STRAIN NAMING above."""
    return LABELS.get(_canonical(name), name)


# ── rcParams ──────────────────────────────────────────────────
# One size ladder, derived from BASE_PT, so changing the base rescales the
# whole figure family coherently. Four rungs, each with a distinct job —
# anything drawn at a size not on this ladder reads as a different figure.
BASE_PT = 8.5

SIZE_TITLE = BASE_PT + 1.0   # figure title, panel letters (both bold)
SIZE_BODY = BASE_PT          # panel titles (bold) and axis labels (regular)
SIZE_TICK = BASE_PT - 1.0    # tick labels, legend entries, footnotes
SIZE_SMALL = BASE_PT - 2.0   # in-figure annotations (per-bar n, and the like)

# One weight for every interval in the repository, a shade under the series
# line: an interval qualifies its marker, it is not a second series.
ERRORBAR_LW = 0.8

RC_PARAMS = {
    # Figure
    "figure.figsize": (COL_DOUBLE, 3.2),
    "figure.facecolor": SURFACE,
    "figure.edgecolor": SURFACE,
    "figure.dpi": 150,
    "figure.titlesize": SIZE_TITLE,
    "figure.titleweight": "bold",
    "figure.labelsize": SIZE_BODY,

    # Type
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
    "font.size": BASE_PT,
    "text.color": INK,
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "mathtext.default": "regular",

    # Axes
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.7,
    "axes.titlesize": SIZE_BODY,
    "axes.titleweight": "bold",
    "axes.titlecolor": INK,
    "axes.titlepad": 5.0,
    "axes.labelsize": SIZE_BODY,
    "axes.labelcolor": INK,
    "axes.labelpad": 3.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.axisbelow": True,
    "axes.prop_cycle": cycler(color=SERIES_CYCLE),

    # Grid — a hairline, behind the data, y only unless a plot asks otherwise
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": GRID,
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
    "grid.alpha": 1.0,

    # Ticks
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "xtick.labelsize": SIZE_TICK,
    "ytick.labelsize": SIZE_TICK,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 2.8,
    "ytick.major.size": 2.8,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "xtick.major.pad": 2.5,
    "ytick.major.pad": 2.5,

    # Marks — the marker carries the datum, the line only joins them, so the
    # line stays well under the marker's width. At column width a line much
    # above 1 pt swallows both the marker and the error bar it sits on.
    "lines.linewidth": 0.9,
    "lines.markersize": 4.2,
    "lines.markeredgewidth": 0.7,
    "lines.solid_capstyle": "round",
    "lines.dash_capstyle": "round",
    "patch.linewidth": 0.6,
    "errorbar.capsize": 2.0,
    "scatter.edgecolors": "none",

    # Legend
    "legend.frameon": False,
    "legend.fontsize": SIZE_TICK,
    "legend.title_fontsize": SIZE_TICK,
    "legend.handlelength": 1.5,
    "legend.handletextpad": 0.5,
    "legend.columnspacing": 1.1,
    "legend.labelspacing": 0.35,
    "legend.borderaxespad": 0.4,
    "legend.borderpad": 0.3,

    # Export — vector text stays text
    "savefig.dpi": 600,
    "savefig.facecolor": SURFACE,
    "savefig.edgecolor": SURFACE,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def use() -> None:
    """Install the house style. Idempotent — call once per driver."""
    plt.rcParams.update(RC_PARAMS)


# ── Axis / figure finishing ───────────────────────────────────

def finish_axis(ax, grid: str = "y") -> None:
    """Apply the house chrome to one axis.

    `grid` is "y" (default), "x", "both" or "none". Everything else here is
    already in the rcParams; this re-applies it so per-axis overrides made
    while drawing do not persist.
    """
    ax.set_axisbelow(True)
    if grid == "none":
        ax.grid(False)
    else:
        ax.grid(False)
        ax.grid(True, axis=grid, color=GRID, linestyle="-", linewidth=0.5,
                alpha=1.0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.7)
        ax.spines[side].set_color(AXIS)


def fit_xticklabels(ax, rotation: float = 30.0, pad_px: float = 2.0) -> bool:
    """Rotate the x tick labels, but only if they would otherwise collide.

    How many labels fit depends on the panel count and the label text, so
    measure the laid-out labels and rotate only on actual overlap. Returns
    True if the labels were rotated.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer)
             for t in ax.get_xticklabels() if t.get_text()]
    if len(boxes) < 2:
        return False
    boxes.sort(key=lambda b: b.x0)
    if not any(boxes[i].x1 + pad_px > boxes[i + 1].x0
               for i in range(len(boxes) - 1)):
        return False
    for text in ax.get_xticklabels():
        text.set_rotation(rotation)
        text.set_ha("right")
        text.set_rotation_mode("anchor")
    return True


def legend(ax, *args, **kwargs):
    """Legend with the house defaults (frameless, tight, secondary ink)."""
    kwargs.setdefault("frameon", False)
    leg = ax.legend(*args, **kwargs)
    if leg is not None:
        for text in leg.get_texts():
            text.set_color(INK_SECONDARY)
        if leg.get_title() is not None:
            leg.get_title().set_color(INK)
    return leg


def legend_outside(ax, *args, x: float = 1.02, y: float = 1.0, **kwargs):
    """House legend in its own column to the right of `ax`.

    The key never covers data or a spine. Both layout engines measure an
    axes legend, so the column comes out of the drawing area rather than out
    of the figure's printed width — give the figure the extra width when the
    panel would otherwise be squeezed. On a multi-panel figure, anchor to the
    rightmost axis so the legend clears every panel.
    """
    kwargs.setdefault("loc", "upper left")
    kwargs.setdefault("bbox_to_anchor", (x, y))
    kwargs.setdefault("borderaxespad", 0.0)
    return legend(ax, *args, **kwargs)


def footnote(fig, text: str, y: float = -0.02, size: float = None) -> None:
    """Italic methods note under the figure, in secondary ink.

    Wrapped to the figure's own width, so a long note grows downward into the
    tight bounding box instead of off the sides of the page.
    """
    # ~19 characters per inch at SIZE_TICK italic in the house sans.
    width = max(40, int(fig.get_figwidth() * 19))
    fig.text(0.5, y, textwrap.fill(text.strip(), width=width),
             ha="center", va="top", linespacing=1.4,
             fontsize=size or SIZE_TICK, fontstyle="italic",
             color=INK_SECONDARY)


def reference_line(ax, value: float, *, axis: str = "y", label_text: str = None,
                   where: str = "left"):
    """A neutral rule at a meaningful value of the plotted quantity — a 1:1
    ratio, a zero difference.

    The rule is drawn behind the data in the muted mark colour; its label,
    like every other annotation, is secondary ink. `label_text` is available
    but unused by default: assay limits (detection floors and the like)
    belong in the methods text, not on the axes.
    """
    draw = ax.axhline if axis == "y" else ax.axvline
    line = draw(value, color=INK_MUTED, linestyle=(0, (4, 2)), linewidth=0.8,
                zorder=0)
    if label_text:
        lo, hi = ax.get_xlim() if axis == "y" else ax.get_ylim()
        pos = lo if where == "left" else hi
        ha = "left" if where == "left" else "right"
        if axis == "y":
            ax.text(pos, value, f" {label_text} ", color=INK_SECONDARY, ha=ha,
                    va="bottom", fontsize=SIZE_TICK)
        else:
            ax.text(value, pos, f" {label_text} ", color=INK_SECONDARY, ha=ha,
                    va="bottom", fontsize=SIZE_TICK, rotation=90)
    return line


# ── Export ────────────────────────────────────────────────────
# PNG for screen and drafts; PDF as the vector master for the journal. One
# `save()` call writes both.
SAVE_FORMATS = ("png", "pdf")


_IMAGE_SUFFIXES = (".png", ".pdf", ".svg", ".eps", ".jpg", ".jpeg", ".tif", ".tiff")


def save(fig, path, formats=None, close: bool = True) -> list:
    """Write `fig` to `path` in every configured format. Returns the paths.

    `path` may carry an image extension or none; only the stem is used. The
    suffix is stripped by name rather than with `Path.with_suffix`, so a stem
    containing a dot (a date, a dilution factor) is not truncated.
    """
    path = Path(path)
    name = path.name
    for suffix in _IMAGE_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    path.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in (formats or SAVE_FORMATS):
        out = path.parent / f"{name}.{ext}"
        fig.savefig(out, format=ext)
        written.append(out)
    if close:
        plt.close(fig)
    return written
