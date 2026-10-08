#!/usr/bin/env python3
"""Driver for the supplementary figure covering HS-3's four other σ factors.

    .venv/bin/python scripts/alt_factors_figure.py            # the composed figure
    .venv/bin/python scripts/alt_factors_figure.py --panels   # each panel on its own

Genomes come from the shared assembly cache the rpoS figure fills. The domain
panels need one InterProScan run over the comparators' σ factor proteins, cached
under `data/rpos_sequence/raw_hmmer/sigma_factors__iprscan.json`, so a second
run makes no network call.
"""

from __future__ import annotations

import argparse

import viz as fs

from figures import alt_factors_locus
from seq.rpoS_locus import domains, fetch, figure, io, orthology

# Below this, a full-length alignment is no longer evidence of orthology and
# the panel would be labelling a match the sequence does not make. The
# conventional twilight-zone floor; every pair the four loci draw clears it
# with better than 90% coverage.
IDENTITY_FLOOR = 20.0


def report(cfg, paths, key, arches) -> list:
    """Print what the two panels for one locus rest on. Returns its tracks."""
    print(f"{key}:")
    tracks = io.synteny_tracks(cfg, paths, key)
    for track in tracks:
        anchor, _ = io.locus_protein(cfg, paths, track.tag, key)
        order = "-".join(("ψ" if gene.pseudo else "") +
                         (gene.local or gene.ortholog)
                         for gene in track.genes if gene.ortholog)
        print(f"  {track.label:14s} {anchor.protein_id:18s} {order}")
    for line in domains.summary(cfg, paths, key, arches):
        print(f"    {line}")
    return tracks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", action="store_true",
                        help="also write each panel as its own figure")
    parser.add_argument("--no-fetch", action="store_true",
                        help="fail rather than download a genome or submit a scan")
    args = parser.parse_args()

    fs.use()
    cfg, paths = io.load()
    paths.ensure()
    paths.alt_output_root.mkdir(parents=True, exist_ok=True)
    if not args.no_fetch:
        fetch.ensure_genomes(cfg, paths)

    arches = domains.ensure_scan(cfg, paths, offline=args.no_fetch)

    rows = []
    for key in cfg["alt_factors"]["order"]:
        tracks = report(cfg, paths, key, arches)

        # The core sets lean on RefSeq product strings wherever an annotation
        # gives no gene symbol — σ54 in three genomes, HS-3's ptsN. Aligning
        # each one to its HS-3 counterpart is what turns that match into an
        # orthology call.
        identities = orthology.core_identities(cfg, paths, key)
        rows += [dict(locus=key, **row) for row in identities]
        aligned = {(row["gene"], row["genome"]) for row in identities}
        for track in tracks:
            if track.tag == "hs3":
                continue
            for gene in track.genes:
                if gene.ortholog and (gene.ortholog, track.label) not in aligned:
                    why = ("pseudogene, drawn open and marked ψ"
                           if gene.pseudo else "no protein in the annotation")
                    print(f"  {gene.ortholog} / {track.label}: "
                          f"not aligned ({why})")
        weakest = min(identities, key=lambda r: r["identity_pct"])
        print(f"  weakest core identity vs HS-3: {weakest['gene']} / "
              f"{weakest['genome']} at {weakest['identity_pct']:.0f}% over "
              f"{weakest['coverage_pct']:.0f}% coverage")
        if weakest["identity_pct"] < IDENTITY_FLOOR:
            raise SystemExit(
                f"{key}: {weakest['gene']} in {weakest['genome']} is only "
                f"{weakest['identity_pct']:.0f}% identical to HS-3's, below "
                f"the {IDENTITY_FLOOR:.0f}% floor this panel's labels rest on")

    written = orthology.write_table(
        rows, paths.alt_output_root / "core_identity.csv")
    print(f"wrote {written}")

    if args.panels:
        figure.render_alt(arches=arches)

    out = alt_factors_locus.plot(
        paths.alt_output_root / "alt_factors_locus.png", arches=arches)
    print("wrote " + ", ".join(str(p) for p in out))


if __name__ == "__main__":
    main()
