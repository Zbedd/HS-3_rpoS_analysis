# rpoS_locus

Genome context for HS-3 *rpoS*, and the structure of the *rpoS*⁻ mutant as
sequenced. Supplies panels A, B and D of the manuscript figure
`figures.rpoS_locus_and_mutant`; panel C uses the shared domain renderer. The same readers serve the
supplement `figures.alt_factors_locus`, which gives the four other σ factors
the two views panels B and C give *rpoS*.

## Run

Use the frozen manuscript inputs from a fresh checkout:

```sh
.venv/bin/python scripts/reproduce.py --only rpos alt-factors
```

For work with the local sequence cache, the individual drivers are:

```sh
.venv/bin/python scripts/rpoS_figure.py            # the composed figure
.venv/bin/python scripts/rpoS_figure.py --panels   # each panel on its own

.venv/bin/python scripts/alt_factors_figure.py           # the supplement, A-H
.venv/bin/python scripts/alt_factors_figure.py --panels  # each panel on its own
```

The genomes are HS-3 plus E. coli, P. aeruginosa and V. cholerae. Downloads
use `seq.assemblies` and a cache shared with the other sequence analyses.
The remaining inputs must already be available: the mutant assembly and the vector are tracked in
`data/rpoS_locus/`, and the domain architectures come from RQ2's
InterProScan cache.

The supplement needs one InterProScan submission of its own, covering the
comparators' σ factor proteins, cached at
`data/rpos_sequence/raw_hmmer/sigma_factors__iprscan.json`. Once that exists
both figures run offline; `--no-fetch` refuses to download or submit anything.

## Modules

| module | holds |
|---|---|
| `config.yaml` | source of truth — accessions, comparator set, σ-factor products, ortholog groups, anchor lengths |
| `paths.py` | every path, resolved from `paths.REPO_ROOT` |
| `io.py` | FASTA / GFF3 / GenBank readers; σ factors, neighbourhood windows, synteny tracks. A locus is a config key, so *rpoS* and the four alternates reach one set of readers |
| `junction.py` | the exact-matching walk that locates the integration, and the layout panel D draws |
| `orthology.py` | pairwise identity behind panel B's ortholog labels, and the supplement's |
| `domains.py` | the supplement's InterProScan pass and its architecture rows |
| `panels.py` | the four panel functions |
| `figure.py` | each panel as its own figure, and the two the supplement composes |

## Outputs

`outputs/rpoS_locus/` — `rpoS_locus_and_mutant_ABCD.{png,pdf}` and
`core_identity.csv`. `--panels` also writes chromosome, synteny, integration
and construct-map figures.

`outputs/alt_factors_locus/` — `alt_factors_locus.{png,pdf}` (panels A–H),
`core_identity.csv` covering all four loci, and eight standalone panels
when `--panels` is given.

See `methods.md` for the citation-ready description and for what the junction
call establishes about the mutant.
