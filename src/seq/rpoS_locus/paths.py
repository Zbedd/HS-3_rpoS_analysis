"""Resolve every filesystem path this workstream uses from a single workdir.

The assembly cache is RQ1's, so a genome fetched here is the same bytes RQ1
and RQ2 read. Only the mutant assembly and the vector are local to this
workstream, and those are tracked inputs rather than a regenerable cache.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import paths as hs3_paths


REPO_ROOT = hs3_paths.REPO_ROOT
DEFAULT_WORKDIR = REPO_ROOT


@dataclass(frozen=True)
class Paths:
    workdir: Path = DEFAULT_WORKDIR

    @property
    def data_root(self) -> Path:
        return self.workdir / "data" / "rpoS_locus"

    @property
    def output_root(self) -> Path:
        return self.workdir / "outputs" / "rpoS_locus"

    @property
    def alt_output_root(self) -> Path:
        """The supplement covering the other four σ factors. A separate folder
        because it is a separate figure, drawn from the same genomes."""
        return self.workdir / "outputs" / "alt_factors_locus"

    @property
    def assemblies_root(self) -> Path:
        """Shared with RQ1 and RQ2 — one copy of each genome on disk."""
        return self.workdir / "data" / "rpos_sequence" / "assemblies"

    @property
    def raw_hmmer(self) -> Path:
        """RQ1's InterProScan cache, keyed on a tag. The supplement's σ factor
        scan lands beside RQ1's own, one directory of InterProScan JSON."""
        return self.workdir / "data" / "rpos_sequence" / "raw_hmmer"

    def assembly_dir(self, tag: str) -> Path:
        return self.assemblies_root / tag

    def genome_fna(self, tag: str) -> Path:
        return self.assembly_dir(tag) / "genomic.fna"

    def genome_gff(self, tag: str) -> Path:
        return self.assembly_dir(tag) / "genomic.gff"

    def data(self, name: str) -> Path:
        return self.data_root / name

    def ensure(self) -> None:
        for p in (self.data_root, self.output_root, self.assemblies_root,
                  self.raw_hmmer):
            p.mkdir(parents=True, exist_ok=True)


def load_config() -> dict:
    import yaml
    with (Path(__file__).parent / "config.yaml").open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)
