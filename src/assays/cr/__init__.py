"""Congo Red assay package.

Modules
-------
config        Constants, sample/media maps, paths, plot toggles.
io            DNG decode, PNG save, JSON cache, CSV r/w.
geometry      Red-ring + dish detection, ring-edge snap.
calibration   Per-angle ring sampling + 2-point calibration.
detection     Colony detection (default + alternatives + overrides).
quantify      SpotMetrics + quantify_spots (OD-chromatic metric).
annotate      Annotation overlay + save_failure_diagnostic.
stats         Hodges-Lehmann estimate, exact Wilcoxon CI, dish pairing.
panels        Panel bodies, each drawing into an axes the caller owns.
figure        The standalone bar and lollipop figures.

`scripts/cr_pipeline.py` is the single driver that composes these
modules. The two-panel levels-over-effect figure lives in
`figures.cr_rpoS_effect`, which draws it from the same `panels`.
"""
