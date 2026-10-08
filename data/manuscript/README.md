# Manuscript inputs

This directory holds the input snapshot used by `scripts/reproduce.py`.
`manifest.json` records file sizes, SHA-256 checksums, provenance, and the
working path expected by each analysis. Existing laboratory inputs in
`data/SURV/` and `data/rpoS_locus/` are included in the manifest at their original
locations.

## Contents

| Files | Source |
|---|---|
| `CR/quantification.csv`, `CR/cr_summary.csv` | Saved colony measurements from the CR pipeline |
| `CR/cache/*.json` | Saved plate geometry, calibration, colony masks, and exclusions |
| `sequence/assemblies/` | RefSeq FASTA, GFF, and protein files for HS-3, E. coli K-12, P. aeruginosa PA14, and V. cholerae |
| `sequence/queries/`, `sequence/hs3_queries/` | The six reference proteins and two HS-3 candidates used in the RpoS/RpoD comparison |
| `sequence/alignments/` | Saved Clustal Omega alignment responses |
| `sequence/interpro/` | Saved InterProScan records for the manuscript proteins |

The CR tables and masks retain all 72 images: both readings and all recorded
media. The manuscript figure uses 72 h, LB and R2A. The same physical plates
were photographed at 48 h and 72 h, but image indices do not establish their
correspondence between readings.

The InterProScan files contain the selected response records from the working
cache, with each retained record unchanged. `rpos_refs__iprscan.json` contains
the eight RpoS/RpoD proteins; `sigma_factors__iprscan.json` contains four other
sigma factors across the four genomes. The original working files and their
hashes are recorded in the manifest. These files are frozen responses from
InterProScan 5.78-109.0, Pfam 38.2, and NCBIfam 19.0. The original service
retrieval dates were not recorded; the snapshot date is recorded separately.

## Raw images

The 24 DNGs required by the manuscript workflow are
`72h_LB_1.DNG` through `72h_LB_12.DNG` and
`72h_R2A_1.DNG` through `72h_R2A_12.DNG`. They total approximately 1.11 GiB
and are kept outside Git. `manifest.json` lists their individual checksums.

Pass the image directory with `--raw-images DIR`; it defaults to `data/CR/`.
The reproduction command links the images into its working directory, or
copies them where symlinks are unavailable. The figure uses raw DNGs for the
colony crops, so the large normalized and annotated PNGs are unnecessary.

## Working copies

```sh
python scripts/reproduce.py --check-inputs
python scripts/reproduce.py --prepare-only
```

Both commands check the required file checksums. The second also prepares a
working copy without calculating results. Use `--only` to select analyses and
avoid requiring inputs for the others.

The default run directory is `outputs/reproduction/`. Its `data/` and
`outputs/CR/` hold working copies of the inputs at the paths the existing
drivers read. All generated figures, tables, and presentation alignments stay
in that run directory. Runs can reuse a directory prepared from the same
snapshot; a different snapshot requires a new `--workdir`.

The reproduction command makes no service requests. It fails on missing or
changed inputs before launching an analysis. PDF timestamps use the snapshot
date through `SOURCE_DATE_EPOCH`; the run manifest records the actual run time.

For an intentional input update, use the individual analysis drivers and
their working caches, review the new inputs, then replace the corresponding
snapshot files and update their manifest entries. Changing a configuration or
the dependency pins also requires updating its recorded checksum. The shared
working cache at `data/rpos_sequence/` remains separate from this snapshot.
