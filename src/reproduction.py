"""Stage frozen manuscript inputs and run the existing drivers offline."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import runpy
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import paths


DRIVERS = {
    "rpos": ("rpoS_figure.py", ["--no-fetch"]),
    "alignments": ("rpoS_alignments.py", []),
    "alt-factors": ("alt_factors_figure.py", ["--no-fetch"]),
    "survival": ("surv_pipeline.py", []),
    "cr": ("cr_pipeline.py", ["--plots-only"]),
}
MANIFEST = paths.SOURCE_ROOT / "data" / "manuscript" / "manifest.json"
MARKER = ".hs3-reproduction.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Not a relative input path: {value!r}")
    return path


def check_inputs(manifest: dict, selected: list, raw_dir: Path) -> list:
    """Check file integrity without loading analyses or computing results."""
    selected = set(selected)
    entries = [f for f in manifest["files"]
               if selected.intersection(f["analyses"])]
    checks = [(paths.SOURCE_ROOT / relative_path(f["path"]), f) for f in entries]
    if "cr" in selected:
        checks += [(raw_dir / relative_path(f["name"]), f)
                   for f in manifest["raw_images"]]
    errors = []
    for path, record in checks:
        if not path.is_file():
            errors.append(f"Missing: {path}")
        elif path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            errors.append(f"Checksum mismatch: {path}")
    if errors:
        if "cr" in selected:
            errors.append("Supply the 24 manuscript DNGs with --raw-images DIR; "
                          "see data/manuscript/README.md.")
        raise ValueError("\n".join(errors))
    return entries


def prepare(manifest: dict, entries: list, selected: list,
            raw_dir: Path, workdir: Path, snapshot_hash: str) -> None:
    """Copy mutable working inputs; leave the tracked snapshot untouched."""
    source_data = paths.SOURCE_ROOT / "data"
    if workdir == paths.SOURCE_ROOT or workdir == source_data or source_data in workdir.parents:
        raise ValueError("Choose a workdir outside the repository's input directories.")
    if raw_dir == workdir / "data" / "CR":
        raise ValueError("Supply raw images from outside the working input copy.")
    marker = workdir / MARKER
    if workdir.exists() and any(workdir.iterdir()):
        if not marker.is_file():
            raise ValueError(f"{workdir} is not an empty or managed reproduction directory.")
        saved = json.loads(marker.read_text())
        if saved.get("snapshot_sha256") != snapshot_hash:
            raise ValueError("This workdir belongs to another input snapshot; "
                             "choose a new --workdir.")

    destinations = {relative_path(f["destination"]) for f in manifest["files"]
                    if f.get("destination")}
    destinations.update(Path("data/CR") / f["name"] for f in manifest["raw_images"])
    destinations.add(Path("data/rpos_sequence/rq2/alignments/sigma70.presentation.fasta"))
    if (workdir / "data").exists():
        unexpected = [p.relative_to(workdir) for p in (workdir / "data").rglob("*")
                      if p.is_file() and p.relative_to(workdir) not in destinations]
        if unexpected:
            raise ValueError("Unexpected files in the working input copy; "
                             "choose a new --workdir: " + ", ".join(map(str, unexpected)))

    workdir.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"snapshot_sha256": snapshot_hash}, indent=2) + "\n")
    for entry in entries:
        if not entry.get("destination"):
            continue
        source = paths.SOURCE_ROOT / relative_path(entry["path"])
        target = workdir / relative_path(entry["destination"])
        if target.is_symlink():
            raise ValueError(f"Working input must not be a symlink: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    if "cr" in selected:
        for entry in manifest["raw_images"]:
            target = workdir / "data" / "CR" / relative_path(entry["name"])
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.is_symlink():
                target.unlink()
            elif target.exists():
                if sha256(target) == entry["sha256"]:
                    continue
                raise ValueError(f"Refusing to replace a changed working raw image: {target}")
            try:
                target.symlink_to(raw_dir / entry["name"])
            except OSError:
                shutil.copyfile(raw_dir / entry["name"], target)


def package_versions() -> dict:
    versions = {}
    for line in (paths.SOURCE_ROOT / "requirements-reproduce.txt").read_text().splitlines():
        if "==" not in line or line.startswith("#"):
            continue
        name = line.split("==", 1)[0]
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def run_driver(name: str) -> None:
    """Execute one driver with Python network access disabled."""
    def offline(event, args):
        if event in {"socket.connect", "socket.getaddrinfo", "socket.sendto"}:
            raise RuntimeError("Manuscript reproduction is offline. "
                               "A required input is missing from the working copy.")

    sys.addaudithook(offline)
    filename, arguments = DRIVERS[name]
    script = paths.SOURCE_ROOT / "scripts" / filename
    sys.argv = [str(script), *arguments]
    runpy.run_path(str(script), run_name="__main__")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", choices=DRIVERS,
                        help="run only the named analyses (default: all five)")
    parser.add_argument("--raw-images", type=Path,
                        default=paths.SOURCE_ROOT / "data" / "CR",
                        help="directory containing the 24 manuscript DNGs")
    parser.add_argument("--workdir", type=Path,
                        default=paths.SOURCE_ROOT / "outputs" / "reproduction",
                        help="isolated run folder (default: outputs/reproduction)")
    parser.add_argument("--check-inputs", action="store_true",
                        help="check required file checksums and exit")
    parser.add_argument("--prepare-only", action="store_true",
                        help="check and stage inputs without running analyses")
    parser.add_argument("--_driver", choices=DRIVERS, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args._driver:
        run_driver(args._driver)
        return

    selected = list(dict.fromkeys(args.only or DRIVERS))
    workdir, raw_dir = args.workdir.resolve(), args.raw_images.resolve()
    try:
        manifest = json.loads(MANIFEST.read_text())
        if manifest["schema_version"] != 1:
            raise ValueError("Unsupported manuscript manifest version.")
        entries = check_inputs(manifest, selected, raw_dir)
        print(f"Inputs checked for: {', '.join(selected)}", flush=True)
        if args.check_inputs:
            return
        snapshot_hash = sha256(MANIFEST)
        prepare(manifest, entries, selected, raw_dir, workdir, snapshot_hash)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"{exc}\n")

    run_record = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot_sha256": snapshot_hash,
        "analyses": selected,
        "python": platform.python_version(),
        "system": platform.system(),
        "machine": platform.machine(),
        "package_versions": package_versions(),
        "source_date_epoch": manifest["source_date_epoch"],
        "completed": [],
        "status": "prepared",
    }
    record_path = workdir / "run_manifest.json"

    def save_record():
        record_path.write_text(json.dumps(run_record, indent=2) + "\n")

    save_record()
    print(f"Working inputs: {workdir}", flush=True)
    if args.prepare_only:
        return

    env = dict(os.environ, HS3_WORKDIR=str(workdir),
               SOURCE_DATE_EPOCH=str(manifest["source_date_epoch"]),
               MPLBACKEND="Agg")
    run_record["status"] = "running"
    save_record()
    for name in selected:
        print(f"Running {name}", flush=True)
        try:
            subprocess.run([sys.executable, "-m", "reproduction", "--_driver", name],
                           env=env, cwd=paths.SOURCE_ROOT, check=True)
        except subprocess.CalledProcessError as exc:
            run_record.update(status="failed", failed_analysis=name)
            save_record()
            parser.exit(exc.returncode, f"{name} failed; see {record_path}\n")
        run_record["completed"].append(name)
        save_record()
    run_record["status"] = "complete"
    save_record()
    print(f"Outputs: {workdir / 'outputs'}", flush=True)


if __name__ == "__main__":
    main()
