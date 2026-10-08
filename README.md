# HS-3 rpoS analysis

Code and analysis inputs for the HS-3 *rpoS* study: Congo red binding,
stress survival, sigma-factor comparisons, and the *rpoS* insertion.

## Running the manuscript analyses

Clone this repository, then use Python 3.9–3.12 and an editable install:

```sh
git clone https://github.com/Zbedd/HS-3_rpoS_analysis.git
cd HS-3_rpoS_analysis
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-reproduce.txt -e .
```

The package versions are frozen in `requirements-reproduce.txt`. The reference
environment used Python 3.9.6 on macOS. Figures prefer Arial; a different
installed font can change text placement.

The Congo red figure needs the 24 raw DNGs for 72 h on LB and R2A. Place these
in `data/CR/`, or give their directory with `--raw-images`:

```sh
python scripts/reproduce.py --raw-images /path/to/CR
```

Raw images are supplied separately. Their filenames and checksums are listed
in [the input manifest](data/manuscript/manifest.json). Everything needed for
the sequence and survival analyses is included:

```sh
python scripts/reproduce.py --only rpos alignments alt-factors survival
```

The command checks the inputs, copies them into `outputs/reproduction/`, and
runs offline. Figures and tables are written under
`outputs/reproduction/outputs/`:

| Analysis (`--only`) | Main output |
|---|---|
| `rpos` | `rpoS_locus/rpoS_locus_and_mutant_ABCD.{png,pdf}` |
| `alignments` | `rpos_sequence/rq2/rpoS/panelB_region_*_msa.png` — four alignment supplements |
| `alt-factors` | `alt_factors_locus/alt_factors_locus.{png,pdf}` |
| `survival` | `SURV/survival_all_treatments.{png,pdf}` and death-rate tables |
| `cr` | `CR/plots/cr_rpoS_effect.{png,pdf}` |

Use `--check-inputs` to check the required files without running an analysis,
or `--workdir DIR` to keep a run elsewhere. The run manifest records the input
snapshot and installed package versions.

## Data and methods

[The input notes](data/manuscript/README.md) describe the saved measurements,
sequence files, and service responses. Congo red reproduction starts from the
saved colony measurements and masks; the image-processing driver is available
separately for reprocessing the DNGs.

Recorded strain names are preserved: `kanR` is the marked control labelled
WT in the figures; recorded `WT` is the unmarked parental isolate.

- [Congo red methods](src/assays/cr/methods.md)
- [Survival data and analysis conventions](data/SURV/README.md)
- [Locus and insertion methods](src/seq/rpoS_locus/methods.md)
- [Sequence alignment methods](src/seq/rpos_alignments/methods.md)

MIC source measurements are not included.

Code is distributed under the [MIT license](LICENSE).

## Publication snapshot

This repository starts with a curated manuscript snapshot and a fresh history.
Contributor credit and the development source revision are recorded in
[PROVENANCE.txt](PROVENANCE.txt).

To export another source snapshot with the manuscript input allowlist:

```sh
python scripts/export_publication.py --destination ../HS3_manuscript
```

The destination must be new and outside this project. This exports the source,
methods, dependency pins, and frozen inputs. Supply the raw CR images separately
as described above. The exported directory can be installed and reproduced on
its own.

## AI disclosure

AI tools assisted with code development, documentation, and validation.
