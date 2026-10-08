"""
Statistics for the SURV stress-survival assay.

Everything is computed on Log10_CFU, the space the figures also plot in.

Families
--------
1. Two-way ANOVA per treatment — Strain, Time, and their interaction.
2. Strain-vs-strain contrasts within each timepoint.
3. Death-rate slopes per strain, compared two ways: pooled ANCOVA over all
   plates, and a summary-measures test over per-replicate slopes.

Multiplicity
------------
Tukey is fitted separately at each timepoint, so `p_tukey` controls the
family-wise error rate across the strain contrasts *at that timepoint* and
nothing else. Because the reported contrasts span all three treatments, the
repository's correction (`multiplicity`) is applied across the whole reported
set — but to the *unadjusted* pairwise p-values (`p_unadj`), never to Tukey's
already-adjusted ones. Each column therefore has one clean interpretation:

    p_tukey                  FWER-controlled within its timepoint
    p_unadj                  no correction
    p_unadj_bonferroni       FWER-controlled across all reported contrasts

Timepoints where every plate is a non-detect are excluded: after imputation
they hold a single repeated constant and support no contrast.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import levene, linregress, t as t_dist, ttest_ind
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import multiplicity

from . import config, io


def significance_label(p):
    """Conventional star notation for a p-value."""
    if pd.isna(p):
        return "na"
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def _strains_present(subset):
    present = set(subset["Strain"].astype(str))
    return [s for s in config.STRAIN_ORDER if s in present]


# ── Two-way ANOVA ─────────────────────────────────────────────

def run_anova(df, treatment_name):
    """Two-way ANOVA of Log10_CFU on Strain, Time and their interaction.

    Time enters as a factor here: the question is whether strains differ at
    the timepoints assayed, not whether the response is linear in time.

    Type II sums of squares, appropriate for the mildly unbalanced design.
    When the interaction is significant the Strain main effect is an
    unweighted average across timepoints and should not be read as "the
    strains differ" — `interaction_significant` flags that case, and the
    within-timepoint contrasts are the interpretable comparison.

    Levene's test across the Strain x Time cells is reported alongside,
    because both the F test and the pooled-error contrasts assume equal
    variances.
    """
    subset = df[df["Treatment"] == treatment_name]
    if subset["Strain"].nunique() < 2:
        return None

    model = ols("Log10_CFU ~ C(Strain) * C(Time)", data=subset).fit()
    table = sm.stats.anova_lm(model, typ=2)

    cells = [g["Log10_CFU"].values for _, g in
             subset.groupby(["Strain", "Time"], observed=True)
             if len(g) > 1]
    levene_p = float(levene(*cells, center="median").pvalue) if len(cells) > 1 else np.nan

    p_interaction = table.loc["C(Strain):C(Time)", "PR(>F)"]
    return {
        "Treatment": treatment_name,
        "p_strain": table.loc["C(Strain)", "PR(>F)"],
        "p_time": table.loc["C(Time)", "PR(>F)"],
        "p_interaction": p_interaction,
        "interaction_significant": bool(p_interaction < config.ALPHA),
        "levene_p_equal_variance": levene_p,
        "mse_resid": model.mse_resid,
        "n_obs": int(model.nobs),
        "df_resid": int(model.df_resid),
    }


# ── Within-timepoint strain contrasts ─────────────────────────

def run_tukey(df, treatment_name):
    """Strain-vs-strain contrasts at each timepoint.

    Tukey is fitted per timepoint rather than across all Strain x Time groups.
    Fitting across all of them would adjust every p-value for the 105
    cross-timepoint comparisons and then discard them, leaving the retained
    contrasts carrying an adjustment for tests that were never reported.

    `mean_diff_log10` is mean(Strain_2) - mean(Strain_1), in log10 CFU/mL.
    Numbers are taken from the result arrays rather than the rendered summary
    table, which rounds to four decimals.
    """
    subset = df[df["Treatment"] == treatment_name].copy()
    if subset["Strain"].nunique() < 2:
        return pd.DataFrame()

    excluded = io.fully_censored_groups(subset)
    skip_times = set(excluded["Time"]) if not excluded.empty else set()

    rows = []
    for time, cell in subset.groupby("Time", observed=True):
        if int(time) in skip_times or cell["Strain"].nunique() < 2:
            continue

        strains = cell["Strain"].astype(str)
        tukey = pairwise_tukeyhsd(endog=cell["Log10_CFU"], groups=strains,
                                  alpha=config.ALPHA)
        labels = pd.DataFrame(tukey.summary().data[1:],
                              columns=tukey.summary().data[0])
        confint = np.asarray(tukey.confint, dtype=float)

        for i in range(len(labels)):
            s1, s2 = str(labels["group1"][i]), str(labels["group2"][i])
            a = cell.loc[strains == s1, "Log10_CFU"]
            b = cell.loc[strains == s2, "Log10_CFU"]
            # Welch, unadjusted — the input to the across-family corrections.
            # Welch rather than pooled because cell variances differ widely.
            p_unadj = float(ttest_ind(b, a, equal_var=False).pvalue)

            rows.append({
                "Treatment": treatment_name,
                "Time": int(time),
                "Strain_1": s1,
                "Strain_2": s2,
                "n_1": int(len(a)),
                "n_2": int(len(b)),
                "mean_diff_log10": float(tukey.meandiffs[i]),
                "ci_lower": confint[i, 0],
                "ci_upper": confint[i, 1],
                "p_tukey": float(tukey.pvalues[i]),
                "p_unadj": p_unadj,
            })

    if not rows:
        return pd.DataFrame()
    return (pd.DataFrame(rows)
            .sort_values(["Time", "Strain_1", "Strain_2"])
            .reset_index(drop=True))


# ── Death-rate slopes ─────────────────────────────────────────

def calculate_death_slopes(df, treatment_name):
    """Log-linear death rate per strain, with a 95% CI on the slope.

    x = time, y = log10(CFU/mL); the slope is the log10 reduction per minute.
    Censored plates are included at their imputed value, which flattens the
    tail where a strain reaches the detection floor; `n_censored` records how
    many points that affects.
    """
    subset = df[df["Treatment"] == treatment_name]
    rows = []

    for strain in _strains_present(subset):
        strain_data = subset[subset["Strain"].astype(str) == strain]
        x = strain_data["Time"].astype(float)
        y = strain_data["Log10_CFU"].astype(float)
        dof = len(x) - 2
        if dof < 1 or x.nunique() < 2:
            continue

        result = linregress(x, y)
        tcrit = t_dist.ppf(1 - config.ALPHA / 2, dof)
        rows.append({
            "Treatment": treatment_name,
            "Strain": strain,
            "slope_log10_per_min": result.slope,
            "slope_ci_lower": result.slope - tcrit * result.stderr,
            "slope_ci_upper": result.slope + tcrit * result.stderr,
            "stderr": result.stderr,
            "r_squared": result.rvalue ** 2,
            "n_obs": int(len(x)),
            "n_censored": int(strain_data["Censored"].sum()),
        })

    return pd.DataFrame(rows)


def replicate_slopes(df, treatment_name):
    """One death-rate slope per biological replicate.

    The unit of replication in this assay is the culture, not the plate.
    Fitting a slope per replicate reduces each culture to a single number,
    which `compare_death_slopes` then compares across strains without
    assuming plates within a culture are independent.
    """
    subset = df[df["Treatment"] == treatment_name]
    rows = []
    for (strain, replicate), g in subset.groupby(["Strain", "Replicate"],
                                                 observed=True):
        if g["Time"].nunique() < 3:
            continue
        result = linregress(g["Time"].astype(float), g["Log10_CFU"].astype(float))
        rows.append({
            "Treatment": treatment_name, "Strain": str(strain),
            "Replicate": int(replicate),
            "slope_log10_per_min": result.slope,
            "r_squared": result.rvalue ** 2,
            "n_obs": int(len(g)),
        })
    return pd.DataFrame(rows)


def compare_death_slopes(df, treatment_name):
    """Pairwise test of whether two strains die at different rates.

    Two tests per pair, because they make different independence assumptions
    and can disagree:

    p_ancova_pooled
        Log10_CFU ~ C(Strain) * Time over every plate; the interaction term
        is the slope difference. Treats plates as independent, so replicate
        structure is ignored and the test can be anticonservative.
    p_replicate_level
        Welch t-test on the per-replicate slopes from `replicate_slopes`.
        The unit of analysis is the culture (n = 4 per strain), which matches
        how the experiment was replicated.

    Where the two disagree the replicate-level test is the conservative
    reading; both are reported rather than one being chosen silently.
    """
    subset = df[df["Treatment"] == treatment_name]
    strains = _strains_present(subset)
    per_rep = replicate_slopes(df, treatment_name)
    rows = []

    for s1, s2 in itertools.combinations(strains, 2):
        pair_data = subset[subset["Strain"].astype(str).isin([s1, s2])].copy()
        pair_data["Strain"] = pair_data["Strain"].astype(str)

        model = ols("Log10_CFU ~ C(Strain) * Time", data=pair_data).fit()
        table = sm.stats.anova_lm(model, typ=2)
        interaction = [i for i in table.index if ":" in i]
        if not interaction:
            continue
        p_pooled = table.loc[interaction[0], "PR(>F)"]

        a = per_rep.loc[per_rep["Strain"] == s1, "slope_log10_per_min"]
        b = per_rep.loc[per_rep["Strain"] == s2, "slope_log10_per_min"]
        if len(a) > 1 and len(b) > 1:
            p_rep = float(ttest_ind(a, b, equal_var=False).pvalue)
            diff = float(b.mean() - a.mean())
        else:
            p_rep, diff = np.nan, np.nan

        rows.append({
            "Treatment": treatment_name,
            "Strain_1": s1,
            "Strain_2": s2,
            "slope_diff_log10_per_min": diff,
            "p_ancova_pooled": p_pooled,
            "p_replicate_level": p_rep,
            "n_replicates_1": int(len(a)),
            "n_replicates_2": int(len(b)),
            "n_plates": int(model.nobs),
        })

    return pd.DataFrame(rows)


# ── Multiplicity ──────────────────────────────────────────────

def apply_corrections(table, p_column, methods=None):
    """Add adjusted p-values across every row of a family.

    Applied over the whole reported set — all treatments together — so the
    correction reflects the number of tests actually reported. Rows with a
    missing p-value are carried through as NaN rather than being allowed to
    void the entire adjusted column.

    The star column is derived from the adjusted value, not the raw one, so a
    reader scanning it sees the corrected call.
    """
    methods = methods or config.CORRECTION_METHODS
    out = table.copy()
    if out.empty:
        return out

    pvals = out[p_column].astype(float).to_numpy()
    for method in methods:
        out[f"{p_column}_{method}"] = multiplicity.correct(pvals, method)

    out[f"{p_column}_raw_label"] = [significance_label(p) for p in pvals]
    out[f"{p_column}_{multiplicity.METHOD}_label"] = [
        significance_label(p)
        for p in out[f"{p_column}_{multiplicity.METHOD}"]]
    return out


# ── Replicate structure ───────────────────────────────────────

def replicate_diagnostic(df):
    """Does replicate identity explain any of the within-cell variation?

    `run_anova` treats plates as independent within a Strain x Time cell.
    This refits with a C(Replicate) block term and reports the share of
    residual *sum of squares* it accounts for — a quantity bounded in [0, 100]
    — alongside its p-value. The change in mean square is also reported; that
    one can go negative, because the block term costs degrees of freedom
    whether or not it explains anything.
    """
    rows = []
    for treatment in [t for t in config.TREATMENT_ORDER if t in set(df["Treatment"])]:
        g = df[df["Treatment"] == treatment]
        base = ols("Log10_CFU ~ C(Strain) * C(Time)", data=g).fit()
        blocked = ols("Log10_CFU ~ C(Strain) * C(Time) + C(Replicate)", data=g).fit()
        table = sm.stats.anova_lm(blocked, typ=2)

        rows.append({
            "Treatment": treatment,
            "p_replicate": table.loc["C(Replicate)", "PR(>F)"],
            "ss_share_of_residual_pct": 100 * (1 - blocked.ssr / base.ssr),
            "mse_change_pct": 100 * (blocked.mse_resid / base.mse_resid - 1),
            "mse_without_block": base.mse_resid,
            "mse_with_block": blocked.mse_resid,
        })
    return pd.DataFrame(rows)


# ── Driver ────────────────────────────────────────────────────

def run_stats(df, treatments=None):
    """Run every family over every treatment and return the assembled tables.

    Returns a dict of DataFrames: anova, tukey, slopes, replicate_slopes,
    slope_comparisons, replicates, excluded_groups. Multiplicity corrections
    are applied across treatments before the tables are returned.
    """
    treatments = treatments or [t for t in config.TREATMENT_ORDER
                                if t in set(df["Treatment"])]

    anova_rows = []
    parts = {"tukey": [], "slopes": [], "replicate_slopes": [],
             "slope_comparisons": [], "excluded_groups": []}

    for treatment in treatments:
        result = run_anova(df, treatment)
        if result:
            anova_rows.append(result)
        parts["tukey"].append(run_tukey(df, treatment))
        parts["slopes"].append(calculate_death_slopes(df, treatment))
        parts["replicate_slopes"].append(replicate_slopes(df, treatment))
        parts["slope_comparisons"].append(compare_death_slopes(df, treatment))
        parts["excluded_groups"].append(
            io.fully_censored_groups(df[df["Treatment"] == treatment]))

    def _concat(frames):
        frames = [f for f in frames if not f.empty]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    tables = {name: _concat(frames) for name, frames in parts.items()}
    tables["anova"] = pd.DataFrame(anova_rows)
    tables["replicates"] = replicate_diagnostic(df)

    if not tables["tukey"].empty:
        tables["tukey"] = apply_corrections(tables["tukey"], "p_unadj")
    if not tables["slope_comparisons"].empty:
        tables["slope_comparisons"] = apply_corrections(
            tables["slope_comparisons"], "p_replicate_level")

    return tables


def lod_sensitivity(load_fn, divisors=None, primary=None):
    """Repeat every test across the LOD imputation conventions.

    Returns one frame per family, tagged by divisor, so a contrast whose call
    depends on where non-detects are placed can be identified directly. The
    divisor used for the primary analysis is always included.
    """
    divisors = list(divisors or config.LOD_SENSITIVITY_DIVISORS)
    if primary is not None and primary not in divisors:
        divisors.append(primary)

    collected = {}
    for divisor in sorted(set(divisors)):
        for name, table in run_stats(load_fn(lod_divisor=divisor)).items():
            if table.empty:
                continue
            tagged = table.copy()
            tagged.insert(0, "lod_divisor", divisor)
            collected.setdefault(name, []).append(tagged)

    return {name: pd.concat(frames, ignore_index=True)
            for name, frames in collected.items()}
