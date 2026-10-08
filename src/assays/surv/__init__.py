"""
SURV — stress-survival (killing curve) assay.

Heat, H2O2 and EtOH exposure of HS-3 WT, kanR and rpoS::ΩkanR — the recorded
names, drawn as Parental, WT and rpoS- — tracked as CFU/mL over time.

    config  paths, strain order, LOD convention, plotting style
    io      load_survival(): CFU/mL, detection limit, censoring
    stats   ANOVA, within-timepoint contrasts, death-rate slopes, multiplicity
    figure  killing-curve figures

`figure` is not imported here: it selects the Agg matplotlib backend at import
time, which should not happen as a side effect of reading the config.

Driven by scripts/surv_pipeline.py.
"""

from __future__ import annotations

from . import config, io, stats

__all__ = ["config", "io", "stats"]
