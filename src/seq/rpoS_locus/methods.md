# Methods — the HS-3 *rpoS* locus and its disruption

Sources and calculations for the A–D manuscript figure. Panel C uses the
RpoS domain analysis described in `seq/rpos_alignments/methods.md`.

Regenerate with:

```sh
.venv/bin/python scripts/reproduce.py --only rpos
```

## Genomes

HS-3 and three reference organisms are represented by the following RefSeq
assemblies. `seq.assemblies` retrieves the sequences and annotations through
the NCBI Datasets API into the shared assembly cache:

| Tag | Assembly | Replicon | Organism |
|---|---|---|---|
| hs3 | GCF_015140455.1 | NZ_AP024094.1 | *Jeongeupia sacculi* HS-3 |
| k12 | GCF_000005845.2 | NC_000913.3 | *Escherichia coli* K-12 MG1655 |
| pa14 | GCF_000014625.1 | NC_008463.1 | *Pseudomonas aeruginosa* UCBPP-PA14 |
| vch | GCF_008369605.1 | NZ_CP043554.1 | *Vibrio cholerae* O1 El Tor N16961 |

The reference set is one genome per bacterial order — Enterobacterales,
Pseudomonadales and Vibrionales, all Gammaproteobacteria; HS-3 is a
Betaproteobacterium of the order Neisseriales. The three are the organisms of
the panel RQ2 uses for the residue-level RpoS call, so the domain panel and the
synteny panel describe the same organisms. Panel C's *V. cholerae* RpoS is
RQ2's UniProt O51804, which differs from the RpoS annotated in N16961
(WP_000116735.1) at one residue.

## Panel A — σ factor complement

Panel A shows five factors: RpoD (WP_200915120.1), RpoS
(WP_255517369.1), RpoH (WP_200916001.1), RpoE (WP_200917341.1) and
RpoN / σ54 (WP_200916434.1). CDSs are selected by RefSeq product strings;
σ54 has no `gene=` attribute in this annotation. Each configured factor must
match exactly one CDS. The selected factors are listed in the configuration.

The chromosome is 3,396,565 bp and circular. The 2,078 bp plasmid pJHS3
(NZ_AP024095.1) carries no σ factor and is not drawn.

## Panel B — *rpoS* neighbourhood

CDSs overlapping the window from 3,300 bp upstream to 1,500 bp downstream of
*rpoS* are drawn in its orientation. K-12 and *V. cholerae* carry it on the
minus strand, so their windows are flipped. All four conserve
**surE – pcm – nlpD – rpoS**, in that order, on one strand.

Orthology is asserted only for those four genes. RefSeq names *surE* and
*rpoS* in all four genomes, but leaves the *nlpD*-like peptidoglycan
DD-metalloendopeptidase unnamed in HS-3 and PA14, and *pcm* unnamed in HS-3,
PA14 and *V. cholerae*; those are matched on the product string, which finds a
candidate but does not establish orthology. Each candidate is therefore
aligned to its HS-3 counterpart — global Needleman-Wunsch, BLOSUM62, gap
open −11, extend −1 — and the identities are
written to `outputs/rpoS_locus/core_identity.csv`. All twelve comparisons
fall between **41.5 % and 56.5 % identity over 91.9–100 % coverage**; the
lowest is *nlpD* against *V. cholerae*, the highest *pcm* against
*P. aeruginosa*. Genes outside those four are drawn as context and are neither
labelled nor claimed.

Pairwise alignments are calculated with Biopython.

## Panel D — the integration, as sequenced

### The reading frame

RefSeq annotates the *rpoS* CDS from 1,070,446 (WP_255517369.1, 316 aa). That
is not the initiation codon: an in-frame `ATG` sits 18 bp upstream at
1,070,428 with no stop between the two, giving a 969 bp CDS of 322 aa whose
first six residues are `MNDQID`. `config.rpoS_start` carries that coordinate
and `io.rpoS_cds` is the only place it is applied.

The correction is what the disruption primers were designed against, and it
decides how the mutant reads — see below. Panels A and B stay on the RefSeq
annotation, where 18 bp is invisible and the comparator genomes carry
annotations of their own.

### The junction call

The mutant assembly (`data/rpoS_locus/`, one closed contig, 101× assembled)
is mapped onto the reference genome (`NZ_AP024094.1`) by exact matching. A
unique *k*-mer anchors each flank, the match is extended base by base to the first mismatch, and the
vector between them is identified against the sequenced pRE118 searched
doubled, since it is circular. Recombination is exact, so there is no indel to
place and no alignment to score. The call is repeated at *k* = 25, 31 and 41
and must agree at all three; a disagreement would mean a repeat or an assembly
artefact swallowed an anchor.

The result:

- The vector is present as one contiguous block of **6,031 bp**, spanning
  pRE118 2,322 → 2,286 across its origin. The 35 bp at pRE118 2,287–2,321 are
  absent — the cloning site the homology arm replaced.
- **497 bp of chromosome are duplicated**, reference 1,070,431–1,070,927,
  one copy either side of the vector. The segment opens at base 4 of the CDS — the
  base immediately after the start codon — and runs to base 500. It occurs
  twice in the mutant and once in the reference.
- The construct is therefore 6,528 bp, and the strain is a single-crossover
  ("loop-in") integrant of pRE118 carrying an internal *rpoS* fragment.

### The arm boundaries match the primers that made it

The duplication is independently confirmed by the primer design
(`HS3/Allelic Exchange Disruption/Gibson Fragments and PCR/Disruption Primers
pRE118.xlsx`). Each *rpoS* primer is a pRE118 overhang plus an HS-3-annealing
half, and the annealing halves place the amplicon exactly where the sequence
puts the duplication:

| | HS-3-annealing sequence | anneals at | observed arm |
|---|---|---|---|
| F | `AACGATCAAATCGATATGCTGG` | 1,070,431 | 1,070,431 |
| R | `TTGATGACATGCACCGG` | 1,070,911–1,070,927 | 1,070,927 |

Designed amplicon 497 bp, observed duplication 497 bp.

**The arm opens one base past the start codon.** The forward primer begins at
1,070,431, immediately after the `ATG` at 1,070,428 — the placement a
disruption needs, since it is what stops the downstream copy initiating. The
primer's own `ATG`, fifteen bases in, is the internal codon at CDS base 19
that RefSeq mistook for the start.

### What the mutant carries

**… pcm – nlpD – *rpoS*′ – pRE118 – *rpoS* – dapD …**

- The **5′ copy** keeps the initiation codon and the native promoter, and runs
  into vector sequence at CDS base 500 — 165 codons of 322. Both σ⁷⁰ regions
  that recognise the −10 and −35 elements lie well beyond residue 165, so it
  cannot yield a functional σS.
- The **3′ copy** carries CDS bases 4 to the stop, so it is missing its
  initiation codon, and its promoter is 6 kb away on the far side of the
  vector. Panel D annotates it *Lacks start codon*. The first in-frame `ATG`
  available to it is the internal one at CDS base 19.

The *nlpD*-like gene upstream is intact: it ends at 1,070,421, ten bases
before the arm begins, and its full 918 bp CDS is present in the mutant.

**Why an annotation pipeline calls the mutant *rpoS* intact.** Run on the
mutant, an annotator finds the 3′ copy, starts it at the same internal `ATG`
RefSeq used, and emits a protein identical to WP_255517369.1 — the two share
an MD5. That agreement is an artefact of both pipelines choosing the internal
codon; the 322-residue product cannot be made from that copy. Reading an
intact *rpoS* call as evidence the knockout failed is reading the wrong thing.
The genotype is `rpoS::pRE118` — not `rpoS::ΩkanR`, which denotes an omega
interposon that is not what the sequence shows.

`scripts/rpoS_figure.py` re-asserts on every run that the arm is duplicated
exactly twice and that it does not cover the start codon, and stops if either
fails.

### Vector annotation

pRE118 as sequenced in this laboratory (`data/rpoS_locus/pRE118.gbk`, 6,066 bp,
annotated by pLannotate): R6Kγ origin, so replication needs a *pir*⁺ host and
the plasmid is suicidal in HS-3; *oriT* for conjugal transfer; *aph(3′)-Ia*
(KanR) for selection; *sacB* for sucrose counter-selection of the second
crossover. Each of these was independently recovered from the mutant's own
annotation at the insertion site, which supplies the vector features in panel D.

## Supplementary figure — the four other σ factors

`figures.alt_factors_locus` gives *rpoD*, *rpoH*, *rpoE* and *rpoN* the two
views panels B and C give *rpoS*, over HS-3 and the same three reference
organisms in the same row order. The factor selection matches panel A.

Each locus is anchored on its σ factor and oriented on it, exactly as panel B
is. The anchor is matched on its ortholog group rather than on a gene symbol:
only K-12 annotates σ54 as *rpoN*, and HS-3, PA14 and *V. cholerae* carry it
as "RNA polymerase factor sigma-54" with no symbol at all. The conserved order each panel shows is

| locus | conserved order, in the anchor's orientation |
|---|---|
| *rpoD* | *rpsU – dnaG – rpoD*, the macromolecular synthesis operon |
| *rpoH* | *ftsY – ftsE – ftsX – rpoH* |
| *rpoE* | *rpoE – rseA – rseB* |
| *rpoN* | *lptB – rpoN – hpf – ptsN* |

and it holds in all four genomes without exception. *P. aeruginosa* names the
*rpoE* operon after alginate — *algU* (also *algT*) *– mucA – mucB* — and the
panel carries those symbols on that row; the genes are the same three.

*V. cholerae* N16961 (GCF_008369605.1) annotates *dnaG* as a frameshifted
**pseudogene**. It is drawn like any other gene in the row — the panel's claim
is about gene order, not about what each gene encodes — but it is marked ψ,
and it is the one core gene with no identity row, since there is no protein to
align.

Orthology is checked the same way panel B's is — global Needleman-Wunsch,
BLOSUM62, gap open −11, extend −1 — with every core gene aligned to its HS-3
counterpart and the result written to
`outputs/alt_factors_locus/core_identity.csv`. Forty-one of the forty-two
comparisons return a protein (the forty-second is the *V. cholerae* *dnaG*
pseudogene), and all forty-one fall between **25.9 % and 71.4 % identity over
91.9–100 % coverage**. The lowest values belong to *rseB*, *rseA*, *ftsX* and
*ptsN* — a periplasmic anti-σ partner, a membrane permease and a PTS subunit,
none of them well conserved at the sequence level, all of them aligned over
their full length. The check that matters most is HS-3's unnamed *ptsN*, which
is 36.1 % identical to K-12 PtsN, 37.9 % to the *P. aeruginosa* protein and
32.4 % to the *V. cholerae* protein, over at least 99.4 % coverage in each.

Domain architectures of the sixteen σ factor proteins come from one
InterProScan submission (EBI InterProScan 5.78-109.0, Pfam-A 38.2, NCBIfam
19.0 — the pins RQ2 records and the same release that produced every other
InterProScan cache here), and every Pfam-A signature is drawn rather than the
σ70 regions alone. Each panel is aligned on σ70 region 2 (PF04542), the one
region every member of the family carries; RpoH would go ragged on region 1.2,
which K-12 does not score.

The architectures are uniform within each factor and diagnostic between them:

- **RpoD** carries region 1.1 (PF03979) and the non-essential region
  (PF04546) in all four — the Group 1 signature, and the reason panel C's
  Gevin call separates RpoD from RpoS.
- **RpoH** carries regions 1.2, 2 and 4 and no region 1.1.
- **RpoE** carries only regions 2 and 4.2 (PF04542 + PF08281), the ECF
  two-domain architecture.
- **RpoN** carries none of the σ70 regions. It is a σ54 protein — activator
  interacting domain (PF00309), core-binding domain (PF04963) and DNA-binding
  domain (PF04552) — in all four genomes, which is why it takes an ordered
  colour family of its own.
