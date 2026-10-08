"""
Data loading and CFU arithmetic for the SURV assay.

`load_survival` is the single place CFU/mL, the detection limit and the
censoring flag are defined; both the statistics and the figures consume the
frame it returns.

Two quantities are kept apart:

    LOD                  the plate's detection floor — what one colony would
                         have represented at the dilution it was read at.
                         A property of the assay.
    Imputed_CFU_per_mL   where a non-detect is placed below that floor
                         (LOD / config.LOD_DIVISOR). An analytical choice.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

# ── Schema ────────────────────────────────────────────────────

REQUIRED_COLUMNS = [
    "Strain", "Treatment", "Conditions", "Time", "Replicate",
    "CFU_Count", "Dilution_Factor", "Volume_plated_uL",
]

NUMERIC_COLUMNS = [
    "Time", "Replicate", "CFU_Count", "Dilution_Factor", "Volume_plated_uL",
]


# ── CFU arithmetic ────────────────────────────────────────────

def calculate_cfu_ml(count, dilution_factor, volume_uL):
    """CFUs/mL = count * 10**Dilution_Factor / (Volume_plated_uL / 1000)

    Dilution_Factor is recorded as a base-10 exponent, so a value of 4 means
    the plate was read from a 10^4 dilution. The input is an exponent, not
    the dilution multiplier.
    """
    return (count * (10 ** dilution_factor)) / (volume_uL / 1000.0)


def calculate_lod(dilution_factor, volume_uL):
    """Detection floor in CFU/mL for a plate read at a given dilution.

    One colony is the smallest reportable result, so the floor carries the
    same dilution term as the count itself.
    """
    return calculate_cfu_ml(1, dilution_factor, volume_uL)


# ── Loading ───────────────────────────────────────────────────

def _require_columns(df, path):
    """Check the schema before anything tries to index into it."""
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"[surv] {path.name} is missing columns: {missing}")


def _validate(df, path):
    """Fail loudly on anything that would otherwise propagate as NaN or a
    silently dropped group."""
    unknown_strain = sorted(set(df["Strain"]) - set(config.STRAIN_ORDER))
    if unknown_strain:
        raise SystemExit(
            f"[surv] {path.name} contains strains not in config.STRAIN_ORDER: "
            f"{unknown_strain}. Add them to the config or correct the data — "
            f"unlisted strains would be dropped from some tables and not others."
        )
    unknown_treatment = sorted(set(df["Treatment"]) - set(config.TREATMENT_ORDER))
    if unknown_treatment:
        raise SystemExit(
            f"[surv] {path.name} contains treatments not in "
            f"config.TREATMENT_ORDER: {unknown_treatment}."
        )

    for col in NUMERIC_COLUMNS:
        bad = df[df[col].isna()]
        if len(bad):
            raise SystemExit(
                f"[surv] {path.name}: {len(bad)} row(s) have a non-numeric "
                f"{col} (first at row {int(bad.index[0]) + 2}). "
                f"Every plate needs one to compute CFU/mL."
            )
    if (df["Volume_plated_uL"] <= 0).any():
        raise SystemExit(f"[surv] {path.name}: Volume_plated_uL must be > 0.")
    if (df["CFU_Count"] < 0).any():
        raise SystemExit(f"[surv] {path.name}: CFU_Count must be >= 0.")


def load_survival(path=None, lod_divisor=None):
    """Load the survival counts and derive CFU/mL, the LOD and censoring.

    Returns a tidy frame with one row per plate and these derived columns:

        LOD                 detection floor for that plate, CFU/mL
        Censored            True where no colonies were counted
        Imputed_CFU_per_mL  value substituted for a non-detect (LOD / divisor)
        CFU_per_mL          measured value, or the imputed value if censored
        Log10_CFU           log10 of the above; the space the models work in
    """
    path = path or config.DATA_FILE
    lod_divisor = config.LOD_DIVISOR if lod_divisor is None else lod_divisor
    if lod_divisor <= 0:
        raise SystemExit(f"[surv] lod_divisor must be > 0 (got {lod_divisor}).")

    df = pd.read_csv(path, encoding="utf-8")
    df.columns = df.columns.str.strip()
    _require_columns(df, path)

    # Drop wholly blank spreadsheet rows before validating the rest.
    df = df.dropna(subset=["Strain", "Treatment"]).copy()
    for col in ("Strain", "Treatment", "Conditions"):
        df[col] = df[col].astype(str).str.strip()
    for col in NUMERIC_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    _validate(df, path)

    df["LOD"] = calculate_lod(df["Dilution_Factor"], df["Volume_plated_uL"])
    df["Censored"] = df["CFU_Count"] == 0
    df["Imputed_CFU_per_mL"] = df["LOD"] / lod_divisor

    measured = calculate_cfu_ml(df["CFU_Count"], df["Dilution_Factor"],
                                df["Volume_plated_uL"])
    df["CFU_per_mL"] = np.where(df["Censored"], df["Imputed_CFU_per_mL"], measured)
    df["Log10_CFU"] = np.log10(df["CFU_per_mL"])

    df["Strain"] = pd.Categorical(df["Strain"], categories=config.STRAIN_ORDER,
                                  ordered=True)
    df = df.sort_values(["Treatment", "Strain", "Replicate", "Time"])
    return df.reset_index(drop=True)


# ── Censoring summary ─────────────────────────────────────────

def censoring_summary(df):
    """Non-detects per treatment / strain / timepoint, for the run manifest."""
    rows = []
    for (treatment, strain, time), g in df.groupby(
            ["Treatment", "Strain", "Time"], observed=True):
        n_censored = int(g["Censored"].sum())
        if not n_censored:
            continue
        censored = g[g["Censored"]]
        rows.append({
            "Treatment": treatment, "Strain": str(strain), "Time": int(time),
            "n_plates": len(g), "n_censored": n_censored,
            "all_censored": n_censored == len(g),
            "LOD_CFU_per_mL": float(censored["LOD"].min()),
            "imputed_CFU_per_mL": float(censored["Imputed_CFU_per_mL"].min()),
            "n_dilutions_among_censored": int(censored["Dilution_Factor"].nunique()),
        })
    return pd.DataFrame(rows)


def fully_censored_groups(df):
    """(Treatment, Time) cells where every plate of every strain is a non-detect.

    Such a cell carries no measurement: after imputation every observation in
    it is the same constant, so it has zero within-cell variance and supports
    no contrast between strains. Reported so it can be excluded rather than
    presented as a null result.
    """
    rows = []
    for (treatment, time), g in df.groupby(["Treatment", "Time"], observed=True):
        if g["Censored"].all():
            rows.append({"Treatment": treatment, "Time": int(time),
                         "n_plates": len(g)})
    return pd.DataFrame(rows)
