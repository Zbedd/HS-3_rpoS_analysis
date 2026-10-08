"""Does the core cassette hold up as orthology, or only as gene symbols?

Panel B labels four genes across four genomes. Two of them — surE and rpoS —
carry a gene symbol in every annotation. The other two do not: the nlpD-like
peptidoglycan DD-metalloendopeptidase is unnamed in HS-3 and PA14, and pcm is
unnamed in HS-3, PA14 and V. cholerae. Matching those on the RefSeq
product string is enough to find a candidate, not enough to call it an
ortholog, so this module aligns each candidate to its HS-3 counterpart and
reports the identity the label rests on.

Global Needleman-Wunsch with BLOSUM62 and BLAST's affine gap costs (-11/-1),
which is the same scoring the RQ1 searches used; the cached BLAST+ binaries
are a win64 build and do not run on this platform.
"""
from __future__ import annotations

import csv
from pathlib import Path

from . import io
from .paths import Paths


def _aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"
    return aligner


def _identity(aligner, one: str, two: str) -> dict:
    alignment = aligner.align(one, two)[0]
    top, bottom = alignment[0], alignment[1]
    matches = sum(a == b for a, b in zip(top, bottom) if a != "-" and b != "-")
    aligned = sum(1 for a, b in zip(top, bottom) if a != "-" and b != "-")
    return {
        "identity": 100.0 * matches / aligned if aligned else 0.0,
        "aligned": aligned,
        "coverage": 100.0 * aligned / min(len(one), len(two)),
    }


def core_identities(cfg: dict, paths: Paths, key: str = None) -> list:
    """One row per comparator × core gene, against the HS-3 protein.

    `key` names a locus — rpoS's own by default, or one of the four the
    supplement draws, whose product-string fallbacks need the same check.
    """
    aligner = _aligner()
    groups = [g["key"] for g in io.locus(cfg, key)["ortholog_groups"]]

    def proteins(tag):
        _, window = io.neighbourhood(cfg, paths, tag, key)
        faa = io.read_fasta(paths.assembly_dir(tag) / "protein.faa")
        found = {}
        for feature in window:
            group = io.ortholog_of(cfg, feature, key)
            if group and group not in found and feature.protein_id in faa:
                found[group] = (feature, faa[feature.protein_id])
        return found

    hs3 = proteins("hs3")
    rows = []
    for tag in cfg["comparators"]:
        if tag == "hs3":
            continue
        other = proteins(tag)
        for group in groups:
            if group not in hs3 or group not in other:
                continue
            feature, sequence = other[group]
            scores = _identity(aligner, hs3[group][1], sequence)
            rows.append({
                "gene": group,
                "genome": cfg["assemblies"][tag]["short"],
                "hs3_protein": hs3[group][0].protein_id,
                "protein": feature.protein_id,
                "named_in_annotation": bool(feature.gene),
                "identity_pct": round(scores["identity"], 1),
                "aligned_aa": scores["aligned"],
                "coverage_pct": round(scores["coverage"], 1),
            })
    return rows


def write_table(rows: list, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path
