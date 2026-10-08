"""Build presentation MSAs from cached Clustal Omega output.

Reorders rows (HS-3 candidate first within each group, then the references in
config order) and rewrites the FASTA ids to publication-style labels via
``labels.label_from_cached_fasta``.

Cached alignment files in ``data/rpos_sequence/rq2/alignments/`` are never
modified — this module writes new ``*.presentation.fasta`` files alongside.
"""
from __future__ import annotations

from pathlib import Path

from .align import parse_aligned_fasta
from .fetch import get_reference
from .labels import Label, label_from_cached_fasta
from .paths import Paths


def _query_path(paths: Paths, entry: dict) -> Path:
    """Resolve a reference's cached FASTA through the shared fetcher.

    Does not rebuild the path from its parts: the fetcher owns where a
    sequence lives, so label rendering and the alignment read the same file.
    """
    return get_reference(paths, entry)


def build_labels(
    paths: Paths, cfg: dict,
) -> tuple[dict[str, dict[str, Label]], dict[str, Label], dict[str, Label]]:
    """Resolve labels for every reference and HS-3 candidate.

    Returns:
        - ``ref_labels[panel][accession] -> Label``
        - ``hs3_labels[role] -> Label``                    (rpoS, rpoD)
        - ``acc_to_label[acc] -> Label``                   (flat lookup table)
    """
    ref_labels: dict[str, dict[str, Label]] = {}
    flat: dict[str, Label] = {}
    for panel, entries in cfg["references"].items():
        ref_labels[panel] = {}
        for entry in entries:
            path = _query_path(paths, entry)
            lab = label_from_cached_fasta(path, role_hint=panel)
            ref_labels[panel][entry["id"]] = lab
            flat[entry["id"]] = lab

    hs3_labels: dict[str, Label] = {}
    for role, acc in cfg["hs3_candidates"].items():
        # HS-3 protein FASTA was cached by fetch.extract_hs3_protein
        path = paths.queries_root / f"hs3_{acc}.faa"
        lab = label_from_cached_fasta(path, role_hint=role, is_hs3=True)
        hs3_labels[role] = lab
        flat[acc] = lab
    return ref_labels, hs3_labels, flat


def build_sigma70_presentation(
    paths: Paths, cfg: dict,
    ref_labels: dict[str, dict[str, Label]],
    hs3_labels: dict[str, Label],
) -> tuple[Path, list[Label]]:
    """Rewrite the σ⁷⁰ MSA in HS-3-first / grouped order with clean ids.

    Order: HS-3 rpoS → RpoS refs (config order) → HS-3 rpoD → RpoD refs.
    Returns the presentation FASTA path and the ordered list of Label objects.
    """
    src = paths.alignments_root / "sigma70.aln.fasta"
    aln = parse_aligned_fasta(src)

    order: list[tuple[str, Label]] = []
    # HS-3 rpoS first
    hs3_rpoS_acc = cfg["hs3_candidates"]["rpoS"]
    order.append((f"{hs3_rpoS_acc}_HS3rpoS", hs3_labels["rpoS"]))
    # RpoS references in config order
    for entry in cfg["references"]["rpoS"]:
        order.append((f"{entry['id']}_RpoS", ref_labels["rpoS"][entry["id"]]))
    # HS-3 rpoD next
    hs3_rpoD_acc = cfg["hs3_candidates"]["rpoD"]
    order.append((f"{hs3_rpoD_acc}_HS3rpoD", hs3_labels["rpoD"]))
    # RpoD references
    for entry in cfg["references"]["rpoD"]:
        order.append((f"{entry['id']}_RpoD", ref_labels["rpoD"][entry["id"]]))

    out = paths.alignments_root / "sigma70.presentation.fasta"
    with out.open("w", encoding="utf-8") as fh:
        for orig_id, lab in order:
            if orig_id not in aln:
                raise KeyError(f"{orig_id!r} missing from {src.name}")
            fh.write(f">{lab.msa_id}\n{aln[orig_id]}\n")
    return out, [lab for _, lab in order]
