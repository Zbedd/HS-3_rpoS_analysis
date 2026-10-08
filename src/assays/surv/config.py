"""
Configuration for the SURV stress-survival assay.

Source of truth for paths, plotting style, and the analysis conventions the
rest of the package reads. No literals in the analysis modules.
"""

from __future__ import annotations

import multiplicity
import paths
import viz as fs

# ── Paths ─────────────────────────────────────────────────────
DATA_FILE = paths.data("SURV", "survival_data.csv")
OUTPUT_DIR = paths.outputs("SURV")

# ── Assay layout ──────────────────────────────────────────────
# Presentation order, in the data's recorded names: the unmarked parent first,
# then the two kanamycin-marked strains so the rpoS mutant sits next to the
# strain the figures call WT.
STRAIN_ORDER = ["WT", "kanR", "rpoS::ΩkanR"]
TREATMENT_ORDER = ["Heat", "H2O2", "EtOH"]

# Strains drawn in the manuscript figures: the marked control (`kanR`) and
# the targeted rpoS insertion mutant. The manuscript describes separate
# insertion sites for these strains. The unmarked parent is omitted from the
# curves; statistics and every CSV use the full STRAIN_ORDER. Set this to
# STRAIN_ORDER to draw all three.
FIGURE_STRAIN_ORDER = ["kanR", "rpoS::ΩkanR"]

# The exposure named in each panel's y-axis label.
TREATMENT_EXPOSURE = {
    "Heat": "43.5 °C heat",
    "H2O2": "0.6 mM H$_2$O$_2$",
    "EtOH": "20% (v/v) EtOH",
}

# This assay draws one organism's two genotypes, so it keys on genotype rather
# than on strain identity: the control and the mutant take the two hues the CR
# figures are published on, and the two assays' main panels read as a pair.
# Fill is not available for the genotype here — an open marker already means a
# timepoint containing a non-detect (see `figure`).
STRAIN_COLOURS = {s: fs.genotype_colour(s) for s in STRAIN_ORDER}
STRAIN_MARKERS = {s: fs.marker(s) for s in STRAIN_ORDER}

# The paper's names, from `viz` so every assay agrees: the data's `kanR` is
# the figures' WT and the data's `WT` is the Parental strain.
STRAIN_LABELS = {s: fs.label(s) for s in STRAIN_ORDER}

# ── Censoring convention ──────────────────────────────────────
# A plate with no colonies is a non-detect, not a true zero.
#
# The detection floor is a property of the plate — what one colony would have
# represented at the dilution it was read at:
#
#     LOD (CFU/mL) = 10 ** Dilution_Factor / (Volume_plated_uL / 1000)
#
# LOD_DIVISOR is the separate analytical choice of where below that floor the
# non-detect is placed: 2 gives LOD/2, 1 gives the floor itself. Kept in
# distinct columns (LOD, Imputed_CFU_per_mL) since only LOD is a property of
# the assay.
LOD_DIVISOR = 2
LOD_SENSITIVITY_DIVISORS = [1, 2]

# ── Statistics ────────────────────────────────────────────────
# The method is the repository's, not this assay's. What is local is the set
# it runs over: the unadjusted pairwise p-values across the whole reported
# set. Tukey's own p-values are already family-wise corrected within a
# timepoint and are reported separately rather than re-adjusted.
ALPHA = multiplicity.ALPHA
CORRECTION_METHODS = [multiplicity.METHOD]

# ── Plotting style ────────────────────────────────────────────
# The house style owns typography, chrome and export settings; only the
# figure geometry specific to this assay lives here.
RC_PARAMS = fs.RC_PARAMS

# One killing-curve panel; the three-treatment comparison spans the page.
# Both widths are drawn-at-print-size. The key rides inside the axes, in the
# corner a descending curve leaves empty, so neither width pays for a legend
# column.
FIG_SIZE_SINGLE = (fs.COL_ONEHALF, 2.7)
FIG_SIZE_PANEL = (fs.COL_DOUBLE, 2.6)
