"""Export the manuscript source and frozen inputs without development history."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1]
INCLUDE = (
    "README.md", "LICENSE", "PROVENANCE.txt", ".gitignore", "pyproject.toml",
    "requirements.txt", "requirements-reproduce.txt", "src", "scripts",
    "data/manuscript", "data/rpoS_locus", "data/SURV",
)
IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info", ".DS_Store")


def export(destination: Path) -> None:
    destination = destination.expanduser().resolve()
    if destination == SOURCE or SOURCE in destination.parents:
        raise ValueError("Choose a destination outside the source project.")
    if destination.exists():
        raise ValueError("The export destination must not already exist.")
    # Validate before creating anything. A release must contain the pinned inputs.
    manifest = json.loads((SOURCE / "data/manuscript/manifest.json").read_text())
    for record in manifest["files"]:
        data = (SOURCE / record["path"]).read_bytes()
        if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(f"Snapshot checksum mismatch: {record['path']}")
    for name in INCLUDE:
        source = SOURCE / name
        if not source.exists():
            raise ValueError(f"Missing export source: {name}")
        if source.is_symlink() or (source.is_dir() and any(p.is_symlink() for p in source.rglob("*"))):
            raise ValueError(f"Export source contains a symlink: {name}")
    destination.mkdir(parents=True)
    for name in INCLUDE:
        source, target = SOURCE / name, destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if source.is_dir():
            shutil.copytree(source, target, ignore=IGNORE)
        else:
            shutil.copy2(source, target)
    print(f"Manuscript snapshot: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", required=True, type=Path)
    args = parser.parse_args()
    try:
        export(args.destination)
    except ValueError as error:
        parser.error(str(error))
