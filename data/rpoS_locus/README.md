# rpoS_locus — primary inputs

Three files that cannot be regenerated from a public accession, so unlike the
assembly cache under `data/rpos_sequence/` they are tracked.

| File | What it is |
|---|---|
| `mutant_rpoS.fna.gz` | Long-read assembly of HS-3 *rpoS*⁻, one closed 3,403,889 bp contig. |
| `mutant_rpoS_annotation.tsv` | The annotation table shipped with that assembly. |
| `pRE118.gbk` | The integration vector as sequenced, 6,066 bp, annotated by pLannotate. |

## Provenance

Both mutant files come from the Plasmidsaurus run archived at
`HS3/Allelic Exchange Disruption/Sequencing/HS3 rpoS check/X8ZWDZ_2_Sample_2_S/`,
copied from `annotation/X8ZWDZ_2_Sample_2_S.{fna,tsv}`. Reported by that run:
233,876 reads, 1,451,516,901 bp, 426× raw and 101× assembled coverage, one
contig, 3,267 genes.

`pRE118.gbk` is `CRM84M_2_pRE118.gbk` from
`HS3/Allelic Exchange Deletion/Plasmid Sequencing/pRE112 and pRE118/`, the
laboratory's own stock of the vector rather than a database reference — the
35 bp its cloning site differs by is what the homology arm replaced.

The full sequencing runs, including reads and QC, stay in those folders; only
what the figures read is copied here.

## Note

The contig is not in the same coordinate frame or orientation as the
reference genome (`NZ_AP024094.1`): the *rpoS* locus reads on the contig's
minus strand, near 1.66 Mb rather than 1.07 Mb. `seq.rpoS_locus.junction` resolves
the mapping; nothing should assume the frames align.
