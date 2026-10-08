"""How this repository corrects for multiple testing.

One method, everywhere: **Bonferroni** across the whole set of comparisons a
figure or a table reports. An assay states which set that is; it does not
choose the method.

    import multiplicity

    adjusted = multiplicity.correct(pvalues)
    call = adjusted < multiplicity.ALPHA

Bonferroni is the conservative choice, which is what the CR sector layout asks
for: `SECTOR_SAMPLE` is a fixed map at a fixed camera angle, so strain is
confounded with plate position and a marginal effect cannot be separated from
a positional one. No correction repairs that, but the strict one does not
add to it.
"""
from __future__ import annotations

import numpy as np

__all__ = ["ALPHA", "METHOD", "correct"]

ALPHA = 0.05
METHOD = "bonferroni"


def correct(pvalues, method: str = METHOD) -> np.ndarray:
    """Bonferroni-adjusted p-values: each raw p times the number of tests.

    Non-finite entries neither count towards that number nor come back with a
    value, so a family holding an unestimable member is corrected for the
    tests it actually made.
    """
    if method != METHOD:
        raise ValueError(f"this repository corrects by {METHOD!r}, "
                         f"not {method!r}")
    p = np.asarray(pvalues, dtype=float)
    usable = np.isfinite(p)
    out = np.full(p.shape, np.nan, dtype=float)
    out[usable] = np.minimum(1.0, p[usable] * int(usable.sum()))
    return out
