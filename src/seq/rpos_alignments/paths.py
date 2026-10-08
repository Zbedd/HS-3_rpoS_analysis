"""Resolve every filesystem path used by the RQ2 pipeline from a single
workdir. Mirrors the RQ1 layout (data/<topic>/, outputs/<topic>/)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import paths as hs3_paths


REPO_ROOT = hs3_paths.REPO_ROOT
DEFAULT_WORKDIR = REPO_ROOT


@dataclass(frozen=True)
class Paths:
    workdir: Path

    @property
    def data_root(self) -> Path:
        # Shared cache root with RQ1; RQ2-only data sits under rq2/
        return self.workdir / "data" / "rpos_sequence" / "rq2"

    @property
    def output_root(self) -> Path:
        return self.workdir / "outputs" / "rpos_sequence" / "rq2"

    @property
    def output_rpoS(self) -> Path:
        return self.output_root / "rpoS"


    @property
    def queries_root(self) -> Path:
        """RQ2-local cache — holds only sequences extracted from the HS-3
        proteome, which have no accession to fetch by."""
        return self.data_root / "queries"

    @property
    def shared_queries_root(self) -> Path:
        """The one cache for fetched reference proteins, shared with RQ1.

        An accession used by both pipelines resolves to a single file here
        rather than to a per-pipeline copy that can drift.
        """
        return self.workdir / "data" / "rpos_sequence" / "queries"

    @property
    def alignments_root(self) -> Path:
        return self.data_root / "alignments"

    @property
    def raw_hmmer(self) -> Path:
        return self.data_root / "raw_hmmer"

    @property
    def rq1_iprscan_root(self) -> Path:
        return self.workdir / "data" / "rpos_sequence" / "raw_hmmer"

    @property
    def rq1_assemblies_root(self) -> Path:
        return self.workdir / "data" / "rpos_sequence" / "assemblies"

    def ensure(self) -> None:
        for p in (
            self.data_root, self.output_root,
            self.output_rpoS,
            self.queries_root, self.shared_queries_root,
            self.alignments_root, self.raw_hmmer,
        ):
            p.mkdir(parents=True, exist_ok=True)


def load_config(workdir: Path) -> dict:
    import yaml
    cfg_path = Path(__file__).parent / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)
