"""Manuscript figures composed across panels.

    cr_rpoS_effect          Congo Red strain levels over the rpoS effect
    rpoS_locus_and_mutant   the rpoS locus, its conservation, and the mutant

A figure here draws no data of its own. It divides the page with `viz.layout`,
loads through each assay's `io`, and fills each cell with a panel function
from that assay's `panels` module — the same function the assay's own
standalone figure uses, so the two cannot disagree.

A height is a panel's plot box. `stack` puts the gaps in from the style guide
and `fit` gives the page back whatever the labels hanging under an axis take,
so a figure module states no inch it had to measure by eye.

The dependency runs one way: `figures` imports `assays`, never the reverse.
A figure that spans assays therefore has a home that neither assay owns.
"""
