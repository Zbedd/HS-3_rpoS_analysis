"""Where the repository keeps things.

The one place that knows the repo layout. Every module that reads an input
or writes a figure resolves it through here rather than counting `..` from
its own file, so moving a module never silently redirects its I/O.
"""

from __future__ import annotations

import os
from pathlib import Path

# src/paths.py -> src -> repo root. An editable install leaves this file in
# the working tree, so `__file__` still points into the repository.
SOURCE_ROOT = Path(__file__).resolve().parents[1]
# The reproduction driver redirects I/O into an isolated copy of the inputs.
REPO_ROOT = Path(os.environ.get("HS3_WORKDIR", SOURCE_ROOT)).resolve()

DATA_ROOT = REPO_ROOT / "data"
OUTPUT_ROOT = REPO_ROOT / "outputs"


def data(assay: str, *parts: str) -> Path:
    """Input path under `data/<assay>/`."""
    return DATA_ROOT.joinpath(assay, *parts)


def outputs(assay: str, *parts: str) -> Path:
    """Output path under `outputs/<assay>/`. Does not create the directory."""
    return OUTPUT_ROOT.joinpath(assay, *parts)
