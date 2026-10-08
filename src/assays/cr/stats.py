"""The estimator every CR figure draws, and the dish pairing it consumes.

Deliberately free of matplotlib: reading the statistics must not select a
rendering backend as a side effect.
"""
from __future__ import annotations

import numpy as np

from .config import PLOT_MEDIA, TIMES


def _signed_rank_counts(n: int) -> list:
    """Exact null distribution of S+, the sum of positively-signed ranks:
    the coefficients of ∏_{r=1..n} (1 + x^r), by knapsack DP. Index s holds
    how many of the 2^n sign patterns give S+ = s."""
    M = n * (n + 1) // 2
    counts = [0] * (M + 1)
    counts[0] = 1
    for r in range(1, n + 1):
        for s in range(M, r - 1, -1):
            counts[s] += counts[s - r]
    return counts


def signed_rank_p(values: np.ndarray) -> float:
    """Exact two-sided Wilcoxon signed-rank p for `values` against zero.

    Exact, not the normal approximation: the whole null is enumerable at
    these sample sizes. Assumes no zero or tied differences, which the
    continuous OD-chromatic metric satisfies.

    p is floored by n at 2/2^n — 0.000488 at n = 12, 0.000977 at n = 11 —
    so a tier below that floor can never be reached and must not be offered.
    """
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    n = v.size
    if n == 0:
        return float('nan')
    order = np.argsort(np.abs(v), kind='mergesort')
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(1, n + 1)
    s_plus = int(ranks[v > 0].sum())
    counts = _signed_rank_counts(n)
    M = n * (n + 1) // 2
    tail = min(s_plus, M - s_plus)
    return min(1.0, 2.0 * sum(counts[:tail + 1]) / (2 ** n))


def _signed_rank_trim(n: int, alpha: float = 0.05) -> int:
    """Number of extreme Walsh averages to drop from each tail to form
    an exact two-sided (1 − `alpha`) Wilcoxon signed-rank CI on `n`
    paired differences. Returns the largest `d` with P(S+ ≤ d) ≤ α/2
    under the null signed-rank distribution; the CI is then
    [W_sorted[d], W_sorted[M-1-d]], M = n(n+1)/2, W = sorted Walsh
    averages. Returns -1 when `n` is too small for the interval to
    exist at this `alpha` (e.g. n ≤ 5 at α = 0.05). Assumes no zero /
    tied differences — safe for the continuous OD-chromatic metric."""
    M = n * (n + 1) // 2
    counts = _signed_rank_counts(n)
    threshold = alpha / 2 * (2 ** n)
    cum, d = 0, -1
    for s in range(M + 1):
        cum += counts[s]
        if cum <= threshold:
            d = s
        else:
            break
    return d

def wilcoxon_hl_ci(values: np.ndarray, alpha: float = 0.05) -> tuple:
    """Hodges–Lehmann estimate (median of the Walsh averages) with an
    exact two-sided (1 − `alpha`) Wilcoxon signed-rank CI. Returns
    (hl, lo, hi); lo/hi are NaN when the sample is too small to form the
    interval. This is the estimator the signed-rank test inverts, so 'CI
    excludes the baseline' matches the signed-rank test at level `alpha`.

    The one estimator every CR figure draws — raw levels as well as
    paired differences. Nothing here computes a location or an interval
    any other way, so a level and an effect are always the same
    statistic. On raw levels the exact
    coverage assumes the values are symmetric about their centre; on
    paired differences that is the signed-rank test's own assumption.
    """
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    n = v.size
    if n == 0:
        return (np.nan, np.nan, np.nan)
    iu = np.triu_indices(n)                       # i ≤ j (incl. diagonal)
    w = np.sort((v[iu[0]] + v[iu[1]]) / 2.0)      # Walsh averages
    hl = float(np.median(w))
    d = _signed_rank_trim(n, alpha)
    if d < 0:
        return (hl, np.nan, np.nan)
    return (hl, float(w[d]), float(w[w.size - 1 - d]))

def dish_pairs(rows: list, time: str, numerator: str,
                denominator: str) -> list:
    """Per-plate (`plate_index`, denominator, numerator) `OD_chrom`
    triples, one list per media (aligned to `PLOT_MEDIA`), sorted by
    plate. Plates are matched by (media, plate_index) and kept only when
    both samples are present."""
    sub = [r for r in rows
           if r['time'] == time and r['sample'] in (numerator, denominator)]
    by_key = {}
    for r in sub:
        v = r.get('OD_chrom')
        if not isinstance(v, float) or not np.isfinite(v):
            continue
        by_key.setdefault((r['media'], r['plate_index']), {})[r['sample']] = v

    per_media = []
    for media in PLOT_MEDIA:
        triples = []
        for (m, idx), pair in by_key.items():
            if m != media:
                continue
            num = pair.get(numerator)
            den = pair.get(denominator)
            if num is None or den is None:
                continue
            triples.append((idx, den, num))
        per_media.append(sorted(triples))
    return per_media

def estimation_cells(rows: list, comparisons: list,
                      conditions: list) -> list:
    """Per comparison, one dict per condition: the dish pairs, their
    differences, the Hodges–Lehmann estimate + Wilcoxon CI, and the plate
    whose difference sits nearest that estimate."""
    panels = []
    for comp in comparisons:
        num, den = comp['numerator'], comp['denominator']
        by_time = {t: dish_pairs(rows, t, num, den) for t in TIMES}
        cells = []
        for media, time in conditions:
            triples = by_time[time][PLOT_MEDIA.index(media)]
            diffs = np.asarray([n - d for _idx, d, n in triples],
                               dtype=np.float64)
            hl, lo, hi = wilcoxon_hl_ci(diffs)
            rep = None
            if diffs.size:
                target = hl if np.isfinite(hl) else float(np.median(diffs))
                rep = triples[int(np.argmin(np.abs(diffs - target)))][0]
            cells.append({'triples': triples, 'diffs': diffs,
                          'hl': hl, 'lo': lo, 'hi': hi, 'rep': rep})
        panels.append(cells)
    return panels
