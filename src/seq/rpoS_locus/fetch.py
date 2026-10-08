"""Fetch the configured genomes into the shared assembly cache."""
from __future__ import annotations

from seq import assemblies
from .paths import Paths


def ensure_genomes(cfg: dict, paths: Paths, tags=None) -> None:
    for tag in (tags or cfg["comparators"]):
        accession = cfg["assemblies"][tag]["accession"]
        assemblies.ensure_assembly(paths, tag, accession)
        assemblies.ensure_gff(paths, tag, accession)
