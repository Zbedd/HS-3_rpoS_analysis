"""Domain architecture of the σ factor itself, for each locus the supplement
draws.

The proteins come from the same annotation the synteny panel is anchored on,
so a row of one panel and a row of the other are the same gene by
construction. InterProScan runs at EBI — the cached BLAST+ and HMMER binaries
are a win64 build — and the JSON lands in RQ1's cache directory under the tag
`sigma_factors`, one scan covering all four factors across the comparators.

The version pins live in RQ2's config, which is where the toolchain that
produced every other InterProScan cache in the repository is recorded.
"""
from __future__ import annotations

import json

from ..iprscan import scan_candidates
from ..protein_domains import (DomainArchitecture,
                                      architecture_from_iprscan_record,
                                      verify_versions)
from ..rpos_alignments.paths import load_config as rq2_config
from . import io
from .paths import Paths

CACHE_TAG = "sigma_factors"

# Filtered against σ70 region 2, which every σ70-family member carries. The
# profile only decides what `scan_candidates` reports back; the cartoon reads
# the cached JSON, so RpoN scoring nothing against it costs nothing.
TARGET_PROFILE = "PF04542"


def factor_proteins(cfg: dict, paths: Paths) -> dict:
    """{accession: sequence} for every σ factor the supplement draws.

    Keyed on accession because that is what InterProScan echoes back; the
    locus and genome a protein came from are recovered by `architecture_rows`
    resolving the same anchors again.
    """
    proteins = {}
    for key in cfg["alt_factors"]["order"]:
        for tag in cfg["comparators"]:
            feature, sequence = io.locus_protein(cfg, paths, tag, key)
            proteins[feature.protein_id] = sequence
    return proteins


def ensure_scan(cfg: dict, paths: Paths, *, offline: bool = False) -> dict:
    """Every architecture in the cache, keyed on accession.

    Submits to EBI only when the cache does not already cover every protein
    the supplement draws; `offline` fails instead, so a figure rerun can be
    held to disk.
    """
    proteins = factor_proteins(cfg, paths)
    cache = paths.raw_hmmer / f"{CACHE_TAG}__iprscan.json"
    if cache.exists() and cache.stat().st_size:
        try:
            data = json.loads(cache.read_text(encoding="utf-8"))
            covered = {(p.get("xref", [{}])[0].get("id") or "").strip()
                       for p in data.get("results", [])}
            missing = set(proteins) - covered
        except ValueError:
            missing = set(proteins)
        if missing:
            print(f"  [domains] cache is stale (missing {len(missing)}); "
                  f"rebuilding")
            cache.unlink()
    if not cache.exists():
        if offline:
            raise SystemExit(
                f"{cache} does not cover the σ factor panel and --no-fetch "
                f"forbids the InterProScan submission that would rebuild it")
        paths.raw_hmmer.mkdir(parents=True, exist_ok=True)
        scan_candidates(proteins, paths.raw_hmmer, CACHE_TAG, TARGET_PROFILE)

    data = json.loads(cache.read_text(encoding="utf-8"))
    verify_versions(data, cache, rq2_config(paths.workdir))
    return {(record.get("xref", [{}])[0].get("id") or "").strip():
            architecture_from_iprscan_record(record)
            for record in data.get("results", [])}


def architecture_rows(cfg: dict, paths: Paths, key: str,
                      arches: dict) -> list:
    """(display label, DomainArchitecture, is_hs3) for one locus, HS-3 first.

    The row order is `comparators`, which is the row order the synteny panel
    beside it uses, so a genome sits on the same line of both panels.
    """
    rows = []
    for tag in cfg["comparators"]:
        feature, _ = io.locus_protein(cfg, paths, tag, key)
        arch = arches.get(feature.protein_id)
        if arch is None:
            raise SystemExit(f"no InterProScan architecture for {tag} {key} "
                             f"({feature.protein_id})")
        rows.append((cfg["assemblies"][tag]["short"], arch, tag == "hs3"))
    return rows


def domain_span(rows: list) -> int:
    """Longest protein among the rows, for a caller sizing an axis."""
    return max(arch.length for _, arch, _ in rows)


def summary(cfg: dict, paths: Paths, key: str, arches: dict) -> list:
    """One reportable line per genome: accession, length, Pfam architecture."""
    lines = []
    for label, arch, _ in architecture_rows(cfg, paths, key, arches):
        order = sorted(arch.pfams_full, key=lambda a: arch.pfams_full[a][0][0])
        lines.append(f"{label}: {arch.protein_id} {arch.length} aa  "
                     + " ".join(order))
    return lines


__all__ = ["CACHE_TAG", "DomainArchitecture", "architecture_rows",
           "domain_span", "ensure_scan", "factor_proteins", "summary"]
