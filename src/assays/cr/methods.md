# Congo red analysis

`assays.cr` measures colony-associated Congo red binding from raw DNG images.
Settings are in `assays/cr/config.py`. The default figure is
`outputs/CR/plots/cr_rpoS_effect.{png,pdf}`; `--panels` also writes the
standalone colony, levels, and effect panels.

For manuscript reproduction from the frozen measurements and masks, use
`python scripts/reproduce.py --only cr --raw-images /path/to/CR`. Its working
inputs and figure are written under `outputs/reproduction/`. The individual
`cr_pipeline.py` driver remains available for image reprocessing.

## Strains, media, and readings

Each plate is a six-sector dish photographed against an on-plate colour
reference card. Sector 0 is blank; sectors 1–5 contain the recorded samples
`WT`, `kanR`, `rpoS`, `PA14`, and `PA14_rpoS`, respectively.

Recorded names are retained in the CSVs. In figures, `kanR` is labelled WT
(the kanamycin-marked control), `rpoS` is labelled *rpoS*⁻ (the insertion
mutant), and the unmarked `WT` isolate is labelled Parental. The parental
isolate is excluded from every figure and statistical comparison because
its R2A plates were contaminated. Its measurements remain in the CSVs for
provenance. The naming translation is defined in `viz/style.py`.

The figures compare Miller LB and R2A. NaCl-free LB measurements remain in
the CSVs but are excluded from the figures. Missing colonies are recorded
in `NA_SPOTS` and excluded from the paired comparisons.

The same set of physical plates was photographed at 48 h and 72 h, not
necessarily in the same order. Image indices identify dishes within a
reading and do not establish correspondence between readings. The manuscript
uses 72 h only. Summaries at 48 h are kept separate: the analysis neither
pools the readings as independent plates nor pairs them across time by index.

## Image processing and measurement

DNG images are decoded with rawpy into 16-bit linear-light RGB, with unity
gamma, camera white balance, and automatic brightness adjustment disabled.
Colonies are segmented from the calibrated image. Measurements use the raw
linear-light image and a local agar annulus spanning 1.30–1.50 times each
colony's effective radius, offset from the detected mask to follow its shape.
Annotated images support inspection of the masks and reference annuli.

For each channel, `core` is the mean intensity within the colony mask and
`ref` is the median intensity in its annulus. Using natural logarithms:

```text
OD_R     = ln(ref_R / core_R)
OD_G     = ln(ref_G / core_G)
OD_chrom = OD_G − OD_R
         = ln(core_R / core_G) − ln(ref_R / ref_G)
```

This is an apparent chromatic optical density measured in reflectance.
Higher values indicate stronger red relative to green signal compared with
the surrounding agar. The local ratios cancel multiplicative channel gains
shared by the colony and its reference. Raw intensities are used because
black-anchored calibration can produce non-positive values, for which the
logarithm is undefined. Channel values are floored at `EPS = 1e-4`.

## Pairing, estimates, and tests

The dish is the unit of replication. Each mutant is compared with its control
on the same dish: `rpoS − kanR` for HS-3 and `PA14_rpoS − PA14` for PA14.
A pair contributes only when both measurements are available. The manuscript
contains 12 paired dishes per condition, except PA14 in R2A, which has 11.

Levels and paired differences are summarised by the Hodges–Lehmann estimate,
the median of the n(n+1)/2 Walsh averages, with an exact Wilcoxon signed-rank
95% confidence interval. The signed-rank interpretation assumes symmetry;
the implementation assumes no zero or tied absolute differences.

Effects are tested against zero using exact two-sided Wilcoxon signed-rank
p-values. Bonferroni correction covers the four comparisons in the manuscript
figure: two organisms in two media at 72 h. Raw p-values are multiplied by
four and capped at one. No between-organism or between-medium contrast is
tested. Confidence intervals are individual 95% intervals without a
multiplicity adjustment, so an interval excluding zero need not receive a
significance mark.

Marks use adjusted p-values: `***` below 0.001, `**` below 0.01, `*` below
0.05, and `n.s.` otherwise. With 11–12 pairs and four tests, the minimum
adjusted p-values are 0.003906 and 0.001953; `***` is not attainable here.

Strain positions were fixed across dishes. Pairing accounts for variation
between dishes but cannot separate strain effects from systematic effects
of sector position.

## Figure display

Each representative dish is selected by the paired difference nearest its
condition's Hodges–Lehmann estimate. The figure shows the colony crops and
per-pixel maps calculated as `ln(ref_G / G(x,y)) − ln(ref_R / R(x,y))`.
These maps are for display; statistics use the colony-level ratios of mean
intensities. Averaging per-pixel logarithms would give a different quantity.
Photographic crops share tone settings, and the maps share a colour scale
symmetric about zero.

Colour identifies organism; open and filled markers identify control and
mutant. Lines connect observations from the same dish. Each effect arrow
starts at the control's location estimate and spans the paired
Hodges–Lehmann difference, which need not equal the difference between the
two marginal location estimates. Its label is `exp(paired difference)`, the
fold change in the chromatic ratio.
