#!/usr/bin/env python3
"""Driver for the rpoS-locus manuscript figure and its standalone panels.

    .venv/bin/python scripts/rpoS_figure.py            # the four-panel manuscript figure
    .venv/bin/python scripts/rpoS_figure.py --panels   # each panel on its own

Writes `rpoS_locus_and_mutant_ABCD.{png,pdf}`.

Genomes are fetched into the shared assembly cache on first run and read
from disk after that. Everything else — the mutant assembly, the vector, the
InterProScan architectures — is already in the repository or in RQ2's cache,
so a second run makes no network call.
"""

from __future__ import annotations

import argparse

import viz as fs

from figures import rpoS_locus_and_mutant
from seq.rpoS_locus import fetch, figure, io, junction, orthology


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", action="store_true",
                        help="also write each panel as its own figure")
    parser.add_argument("--no-fetch", action="store_true",
                        help="fail rather than download a missing genome")
    args = parser.parse_args()

    fs.use()
    cfg, paths = io.load()
    paths.ensure()
    if not args.no_fetch:
        fetch.ensure_genomes(cfg, paths)

    print("Integration:")
    integration = junction.resolve(cfg, paths, verbose=True)
    rpoS = io.rpoS_cds(cfg, paths)
    split = integration.arm_split(rpoS)
    print(f"  rpoS CDS {rpoS.start:,}..{rpoS.end:,} "
          f"({rpoS.length} bp, {split['cds_residues']} aa)")
    print(f"  arm opens at CDS base {split['first_cds_base']} and covers "
          f"{split['coding']:,} bp ({split['residues']} codons); start codon "
          f"{'INCLUDED' if split['covers_start_codon'] else 'excluded'}")
    if split["covers_start_codon"]:
        raise SystemExit(
            "the arm covers the start codon, so the downstream copy could "
            "initiate; panel D's annotation would be wrong")
    copies = junction.copies_of_arm(cfg, paths, integration)
    if copies != 2:
        raise SystemExit(
            f"the homology arm occurs {copies} times in the mutant, not twice; "
            "the single-crossover model panel D draws does not hold")

    rows = orthology.core_identities(cfg, paths)
    written = orthology.write_table(rows, paths.output_root / "core_identity.csv")
    print(f"wrote {written}")
    weakest = min(rows, key=lambda r: r["identity_pct"])
    print(f"  core cassette vs HS-3: identity "
          f"{weakest['identity_pct']:.0f}%-"
          f"{max(r['identity_pct'] for r in rows):.0f}% "
          f"(weakest {weakest['gene']} / {weakest['genome']})")

    if args.panels:
        figure.render_all()

    paths.output_root.mkdir(parents=True, exist_ok=True)
    out = rpoS_locus_and_mutant.plot(
        paths.output_root / "rpoS_locus_and_mutant_ABCD.png")
    print("wrote " + ", ".join(str(p) for p in out))


if __name__ == "__main__":
    main()
