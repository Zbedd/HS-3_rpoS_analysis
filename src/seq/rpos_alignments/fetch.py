"""Sequence retrieval for RQ2.

Reuses RQ1's UniProt + NCBI fetchers. Extracts HS-3 candidate sequences from
the RQ1 proteome cache (no re-fetch).
"""
from __future__ import annotations

from pathlib import Path

from .. import seqfetch
from .paths import Paths


def get_reference(paths: Paths, entry: dict) -> Path:
    """Resolve one reference protein through the shared fetcher.

    Reads and writes RQ1's cache directory rather than a private RQ2 copy, so
    an accession used by both pipelines (P13445, P04949, WP_003113871.1, ...)
    is fetched once and both pipelines read the same bytes. Previously each
    pipeline kept its own copy, which could drift.
    """
    rec = seqfetch.get_protein(entry["id"], entry.get("source"),
                               cache_dir=paths.shared_queries_root)
    return rec.path


def extract_hs3_protein(paths: Paths, accession: str) -> Path:
    """Pull a single record by accession from the cached HS-3 RefSeq proteome.

    The extracted single-accession file is checked first: once it exists, the
    proteome it came from is no longer needed, so a cached RQ2 regenerates
    without the ~25 MB RQ1 assembly cache on disk.
    """
    out = paths.queries_root / f"hs3_{accession}.faa"
    if out.exists() and out.stat().st_size > 0:
        return out

    proteome = paths.rq1_assemblies_root / "hs3" / "protein.faa"
    if not proteome.exists():
        raise SystemExit(
            f"HS-3 proteome not found at {proteome}. "
            "Run RQ1 first to populate the assembly cache."
        )

    keep = False
    found = False
    lines: list[str] = []
    for line in proteome.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            header_acc = line[1:].split()[0]
            if header_acc == accession:
                keep = True
                found = True
                lines.append(line)
            else:
                keep = False
        elif keep:
            lines.append(line)
    if not found:
        raise SystemExit(
            f"Could not find {accession} in HS-3 RefSeq proteome "
            f"({proteome.name})"
        )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def read_fasta_record(path: Path) -> tuple[str, str]:
    """Return (header_id, sequence) for a single-record FASTA file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = ""
    seq_lines: list[str] = []
    for line in lines:
        if line.startswith(">"):
            header = line[1:].split()[0]
        else:
            seq_lines.append(line.strip())
    return header, "".join(seq_lines)
