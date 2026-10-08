# RpoS alignment methods

Reproduce the frozen inputs with `python scripts/reproduce.py --only alignments`.
The local-cache driver is `python scripts/rpoS_alignments.py`.

## Pipeline overview

`python scripts/rpoS_alignments.py` runs the RpoS/RpoD analysis and writes
its four manuscript alignment supplements and diagnostic tables. Additional
RpoS domain cartoons and logos require `--exploratory`.

1. **Candidate retrieval (no re-fetch).** HS-3 candidate protein
   sequences are extracted from the frozen RefSeq proteome.
   The candidates are `WP_255517369.1` (rpoS), `WP_200915120.1`
   (rpoD; within-genome control).
2. **Reference sequence retrieval.** Three RpoS proteins, three RpoD
   proteins from UniProt KB / NCBI
   RefSeq, all pinned in `config.yaml` and resolved through the shared
   fetcher (`seq.seqfetch`), which validates that each accession
   returns a single non-empty record whose identifier matches the one
   requested.

   The σ⁷⁰ panels sample **one organism per bacterial order** and
   mirror each other organism-for-organism, so a diagnostic column
   reflects the RpoS/RpoD split rather than a difference in taxon
   sampling between the two sides:

   | Organism | Order | RpoS | RpoD |
   |---|---|---|---|
   | *E. coli* K-12 | Enterobacterales | P13445 | P00579 |
   | *P. aeruginosa* PA14 | Pseudomonadales | WP_003113871.1 | WP_003085035.1 |
   | *V. cholerae* | Vibrionales | O51804 | Q9KUK1 |

   All three references are Gammaproteobacteria; HS-3 (*Jeongeupia
   sacculi*) is a Betaproteobacterium of the order Neisseriales. The
   *V. cholerae* RpoS, O51804, differs from the RpoS annotated in the
   N16961 genome (WP_000116735.1) at one residue.
3. **Multiple-sequence alignment.** The cached Clustal Omega alignment contains
   HS-3 RpoS, HS-3 RpoD, three RpoS references, and three RpoD references.
   Clustal text and aligned FASTA are included in the frozen inputs.
4. **Domain architecture (Gevin framework).** Pfam + NCBIfam
   annotation via EBI InterProScan REST v5. The Gevin et al. 2024
   essential-domain rule (*BMC Genomics* 25:512) classifies each
   σ⁷⁰ protein as Group 1 (RpoD-class — must carry **PF03979
   Sigma70_r1.1**) or Group 2 (RpoS-class — must lack PF03979 while
   retaining PF00140 r1.2, PF04542 r2, PF04539 r3, and PF04545 r4
   or PF08281 r4.2).
5. **Column-by-column residue diagnostic.** Walking the σ⁷⁰
   alignment within E. coli RpoS (P13445) residue ranges 57–113
   (region 1.2), 114–186 (region 2), 187–253 (region 3), and 254–322
   (region 4), every alignment
   column where (a) all three RpoS references share residue R_S,
   (b) all three RpoD references share residue R_D, and
   (c) R_S ≠ R_D is classified as a *diagnostic column*. The HS-3
   RpoS candidate's residue in that column is then scored as
   `RpoS-like` (matches R_S), `RpoD-like` (matches R_D), or
   `neither`.
6. **Supplementary figure.** Panel A: domain-architecture cartoon
   rendered programmatically with matplotlib from the InterProScan
   coordinates. Panel B: four trimmed MSA slices, one for each conserved region,
   rendered with `pyMSAviz`; diagnostic columns marked.
7. **Reports.** `RQ2_rpoS_findings.md`,
   `domain_architecture.csv`, `residue_diagnostic.csv`,
   `run_manifest.json`, and the figures.

## Why three independent diagnostics

Sigma factor identity is one of the standard cases where homology
metrics alone do not separate sister genes — RpoS (group 2) and RpoD
(group 1) share full domain architecture in every region except r1.1,
and percent identity over the conserved σ⁷⁰ core does not separate
them either (Chiang & Schellhorn 2010, *J Mol Evol* 70:557–571).
The RQ2 design therefore stacks three orthogonal diagnostics:

1. **Modular domain architecture (Gevin 2024 framework)** — a single
   essential-domain test (PF03979 presence/absence) that captures
   the Group 1 / Group 2 distinction.
2. **Curated family HMMs** — `TIGR02394 / PRK05657.1` for RpoS and
   `TIGR02393 / PRK05658.1` for RpoD. These HMMs are themselves
   trained on the same residue patterns the column-based test
   measures, so the agreement of the curated HMMs with the column
   test is a coherent rather than independent signal — included
   because the curated families are what NCBI itself uses to assign
   the gene symbol.
3. **Column-by-column residue diagnostic** at σ⁷⁰ regions 2 and 4.
   This is the residue-level evidence requested by the validation
   methods document (region 2.4 promoter-recognition residues at
   the -10 element, region 4.2 promoter-recognition residues at
   the -35 element).

If all three agree, the identity call is unambiguous.

## Reproducibility

Paths derive from `--workdir`. Frozen sequence inputs, service responses, and
software pins are checked against `data/manuscript/manifest.json` before
reproduction. The reproduction command blocks network access. Working caches
are staged under `data/rpos_sequence/`; outputs are under
`outputs/rpos_sequence/rq2/rpoS/`. Fresh service requests are not part of
frozen manuscript reproduction.

## Tool stack to cite

EBI Clustal Omega REST (Sievers et al. 2011); EBI InterProScan 5 /
HMMER 3 (Paysan-Lafosse et al. 2023; Eddy 2011); Pfam and NCBIfam HMM
libraries; Gevin et al. 2024, *BMC Genomics* 25:512 (modular domain
identity framework for σ⁷⁰ proteins); Chiang & Schellhorn 2010,
*J Mol Evol* 70:557–571 (RpoS-vs-RpoD reciprocal-best-hit method);
`pyMSAviz` for MSA rendering; matplotlib for the domain cartoon.
