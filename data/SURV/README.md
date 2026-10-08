# SURV — stress-survival (killing curve) data

Time-kill of HS-3 under three stresses, tracked as CFU/mL.

| | |
|---|---|
| Strains | `WT`, `kanR`, `rpoS::ΩkanR` — recorded names. The figures call these Parental, WT and *rpoS*⁻; see `src/README.md`. |
| Treatments | Heat, H2O2, EtOH |
| Replicates | 4 biological per strain × treatment, except EtOH/WT t=60 and Heat/kanR t=60 which have 3 |
| Plates | 178 |

## Files

| File | Contents |
|---|---|
| `survival_data.csv` | Analysis input. Plate-level counts, one row per plate. |
| `survival_data.xlsx` | Archived original workbook the CSV was derived from. |

The workbook preserves the recorded plate counts and assay metadata. The
manuscript workflow uses the CSV derived from it to reproduce CFU calculations,
statistics, and figures. Plate photographs are not included in this checkout;
reproduction starts from the recorded counts.

`survival_data.csv` is regenerated from the workbook with:

```python
df = pd.read_excel("data/SURV/survival_data.xlsx")
df.columns = df.columns.str.strip()
df.dropna(subset=["Treatment", "Strain"]).to_csv(
    "data/SURV/survival_data.csv", index=False, encoding="utf-8")
```

## Columns

| Column | Units | Notes |
|---|---|---|
| `Strain` | — | `WT`, `kanR`, `rpoS::ΩkanR`. `kanR` is the strain the paper calls WT; `WT` here is the unmarked parental isolate. |
| `Treatment` | — | `Heat`, `H2O2`, `EtOH` |
| `Conditions` | mixed | Exposure level as recorded. |
| `Time` | minutes | Heat/H2O2: 0–80. EtOH: 0, 10, 20, 40, 60. |
| `Replicate` | — | 1–4 |
| `CFU_Count` | colonies | Raw count on the plate. `0` is a non-detect. |
| `Dilution_Factor` | log10 exponent | `4` means the plate came from a 10⁴ dilution. |
| `Volume_plated_uL` | µL | 100 throughout. |

**`Dilution_Factor` is a base-10 exponent, not a multiplier.** For example,
`4` denotes a 10⁴ dilution. `assays.surv.io` applies this convention.

## Derived quantities

Defined once, in `assays.surv.io`:

```
CFU/mL = CFU_Count * 10**Dilution_Factor / (Volume_plated_uL / 1000)
LOD    = 1         * 10**Dilution_Factor / (Volume_plated_uL / 1000)
Imputed_CFU_per_mL = LOD / LOD_DIVISOR          # non-detects only
```

`LOD` and `Imputed_CFU_per_mL` are deliberately separate columns. The first is
a property of the plate — what one colony would have represented at the
dilution it was read at. The second is an analytical choice about where to
place a non-detect below that floor, and it is the only one that moves when
`--lod-divisor` changes. Neither is drawn on a figure: the detection floor
belongs in the methods text, and a censored timepoint is marked with an open
marker instead.

A plate with no colonies is a non-detect rather than a true zero, so it is
flagged `Censored` and its `CFU_per_mL` is set to `Imputed_CFU_per_mL`.

Non-detects: 13 plates, all Heat, at t=60 (WT, 1 plate) and t=80 (all three
strains, 4 plates each). All were read at `Dilution_Factor = 1`, so the floor
is 100 CFU/mL and the default LOD/2 imputation places them at 50 CFU/mL.

`LOD_DIVISOR` in `assays.surv.config` selects the convention (default 2).
The pipeline reruns every test across both conventions and writes
`sensitivity_*.csv`, so any result that depends on the choice is visible.

Any timepoint at which *every* plate is a non-detect is excluded from the
strain contrasts and listed in `excluded_groups.csv`. After imputation such a
cell holds one repeated constant, so it has no within-cell variance and
supports no comparison between strains; including it would also deflate the
pooled error used by neighbouring contrasts. Heat t=80 is the one such cell.

## Multiple testing

Three p-value columns appear in `within_timepoint_contrasts.csv`, each with a
single interpretation:

| Column | Controls |
|---|---|
| `p_tukey` | FWER across the strain contrasts **at that timepoint** — Tukey is fitted per timepoint, not across all Strain × Time groups |
| `p_unadj` | nothing; Welch t-test between the two strains |
| `p_unadj_bonferroni` | FWER across all reported contrasts, all treatments |

Bonferroni is the repository's correction throughout (`src/multiplicity.py`).
It is applied to the *unadjusted* p-values, never to Tukey's, so no value is
corrected twice. `p_unadj_bonferroni_label` carries the star notation for the
corrected call; `p_unadj_raw_label` is the uncorrected one.
