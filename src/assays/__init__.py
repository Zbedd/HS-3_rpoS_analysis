"""One package per phenotypic assay.

    cr    Congo Red colony morphology (DNG pipeline -> OD_chrom)
    surv  Stress-survival killing curves

Each exposes the same two things: `io.load…()` for a tidy frame, and
`panels` for functions that draw one panel into an axes the caller owns.
That contract is what lets `figures` compose panels across assays.
"""
