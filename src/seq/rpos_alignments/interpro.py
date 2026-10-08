"""Load and validate cached domain annotations for the local reference panels."""
from __future__ import annotations

import json
from pathlib import Path

from ..iprscan import scan_candidates
from .paths import Paths


from ..protein_domains import (
    DomainArchitecture, verify_versions, architecture_from_iprscan_record, PF_R11,
)


def load_cached_iprscan(paths: Paths, cfg: dict | None = None) -> dict[str, dict]:
    """Load every protein record from RQ1's InterProScan JSON cache.

    If ``cfg`` is provided, each cache file is checked against the pins
    in ``cfg['software_versions']`` before its records are accepted.
    """
    out: dict[str, dict] = {}
    for tag in ("hs3", "pa14", "k12"):
        f = paths.rq1_iprscan_root / f"{tag}__iprscan.json"
        if not f.exists():
            continue
        data = json.loads(f.read_text(encoding="utf-8"))
        if cfg is not None:
            verify_versions(data, f, cfg)
        for prot in data.get("results", []):
            pid = (prot.get("xref", [{}])[0].get("id") or "").strip()
            if pid:
                out[pid] = prot
    return out


def scan_proteins(
    paths: Paths, proteins: dict[str, str], tag: str,
    target_profile: str = PF_R11,
    cfg: dict | None = None,
) -> dict[str, dict]:
    """Submit unscanned proteins to InterProScan and return the raw JSON map.

    Cache validation: if an existing cache exists but does not cover every
    requested accession (e.g., the reference panel changed), the stale cache
    is removed and a fresh full scan is run. Otherwise the cache is reused.

    If ``cfg`` is provided, the cache file is checked against the pins in
    ``cfg['software_versions']`` after loading and before returning.
    """
    if not proteins:
        return {}
    cache = paths.raw_hmmer / f"{tag}__iprscan.json"
    if cache.exists():
        try:
            data = json.loads(cache.read_text(encoding="utf-8"))
            cached_ids = {
                (p.get("xref", [{}])[0].get("id") or "").strip()
                for p in data.get("results", [])
            }
            missing = set(proteins) - cached_ids
            if missing:
                print(f"  [interpro] cache '{tag}' is stale "
                      f"(missing {len(missing)}: {sorted(missing)[:3]}…); rebuilding")
                cache.unlink()
        except Exception:
            cache.unlink()
    scan_candidates(proteins, paths.raw_hmmer, tag, target_profile)
    out: dict[str, dict] = {}
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
        if cfg is not None:
            verify_versions(data, cache, cfg)
        for prot in data.get("results", []):
            pid = (prot.get("xref", [{}])[0].get("id") or "").strip()
            if pid:
                out[pid] = prot
    return out


def load_rq2_iprscan(paths: Paths, cfg: dict | None = None) -> dict[str, dict]:
    """Records for the RQ2 reference panels, read from their cache only.

    The pipeline reaches these through `scan_proteins`, which needs the
    sequences in hand and will submit to EBI when the cache does not cover
    them. A figure only ever wants what is already on disk, so it reads the
    file directly and gets nothing rather than a network call.
    """
    out: dict[str, dict] = {}
    for tag in ("rpos_refs",):
        cache = paths.raw_hmmer / f"{tag}__iprscan.json"
        if not cache.exists():
            continue
        data = json.loads(cache.read_text(encoding="utf-8"))
        if cfg is not None:
            verify_versions(data, cache, cfg)
        for prot in data.get("results", []):
            pid = (prot.get("xref", [{}])[0].get("id") or "").strip()
            if pid:
                out[pid] = prot
    return out


def architecture_rows(paths: Paths, cfg: dict, panels: tuple = ("rpoS",),
                      arches: dict = None) -> list:
    """(label, DomainArchitecture, is_hs3) rows for one or more panels.

    HS-3's own candidate leads each panel, then that panel's references in
    config order. Both RQ2's supplementary cartoon and the manuscript figure
    build their rows here, so the two cannot disagree about which proteins are
    shown or how they are labelled.

    Reads only the caches; no network call.
    """
    from . import msa_view

    if arches is None:
        records = load_cached_iprscan(paths, cfg)
        records.update(load_rq2_iprscan(paths, cfg))
        arches = {pid: architecture_from_iprscan_record(rec)
                  for pid, rec in records.items()}
    selected = {**cfg,
                "references": {p: cfg["references"][p] for p in panels},
                "hs3_candidates": {p: cfg["hs3_candidates"][p] for p in panels}}
    ref_labels, hs3_labels, _ = msa_view.build_labels(paths, selected)

    rows = []
    for panel in panels:
        accession = cfg["hs3_candidates"].get(panel)
        arch = arches.get(accession)
        if arch is None:
            raise SystemExit(
                f"Missing InterProScan architecture for HS-3 {panel} "
                f"({accession}).")
        rows.append((hs3_labels[panel].short, arch, True))
    for panel in panels:
        for entry in cfg["references"][panel]:
            arch = arches.get(entry["id"])
            if arch is not None:
                rows.append((ref_labels[panel][entry["id"]].short, arch, False))
    return rows
