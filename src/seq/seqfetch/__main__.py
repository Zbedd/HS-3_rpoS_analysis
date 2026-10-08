"""Audit accessions declared in the sequence reference config.

    python -m seq.seqfetch                # audit the retained RQ2 config
    python -m seq.seqfetch --config rq2   # explicit config selection
    python -m seq.seqfetch --strict       # exit 1 if anything is suspect

Fetches each accession through the shared fetcher and prints the label the
config claims next to the description the database actually returned, so a
mismatch is visible rather than inferred. Nothing is rewritten automatically:
choosing the right replacement accession is a judgement call.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

import paths

from . import check_label, get_protein

# Resolved against this package's own location, not the working directory:
# the configs are siblings of `seqfetch` and travel with it.
_SEQ = Path(__file__).resolve().parents[1]

CONFIGS = {
    "rq2": (_SEQ / "rpos_alignments" / "config.yaml", "references"),
}


def audit(which: str, workdir: Path) -> list:
    cfg_path, key = CONFIGS[which]
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    panels = cfg.get(key) or {}

    results = []
    shown = cfg_path.relative_to(paths.REPO_ROOT)
    print(f"\n{'=' * 78}\n{which.upper()}  ({shown})\n{'=' * 78}")
    for panel, entries in panels.items():
        print(f"\n  panel: {panel}")
        for entry in entries:
            acc, declared = entry["id"], entry.get("name", "")
            try:
                rec = get_protein(acc, entry.get("source"),
                                  cache_dir=workdir / "data/rpos_sequence/queries")
            except Exception as exc:  # noqa: BLE001 - report, don't abort the audit
                print(f"    [ERROR ] {acc:<18} {declared}")
                print(f"             {exc}")
                results.append(("error", panel, acc, declared, str(exc)))
                continue
            chk = check_label(rec, declared)
            mark = {"ok": "  ok   ", "suspect": "SUSPECT", "unverifiable": "  ?    "}[chk.verdict]
            print(f"    [{mark}] {acc:<18} {declared}")
            print(f"             -> {rec.description[:88]}  ({rec.length} aa)")
            if chk.verdict == "suspect":
                print(f"             !! {chk.reason}")
            results.append((chk.verdict, panel, acc, declared, chk.reason))
    return results


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", choices=[*CONFIGS, "all"], default="all")
    ap.add_argument("--workdir", type=Path, default=paths.REPO_ROOT,
                help="root of the sequence cache (default: the repo root)")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any accession is suspect or errored")
    args = ap.parse_args(argv)

    which = list(CONFIGS) if args.config == "all" else [args.config]
    all_results = []
    for w in which:
        all_results.extend(audit(w, args.workdir))

    bad = [r for r in all_results if r[0] in ("suspect", "error")]
    print(f"\n{'=' * 78}")
    print(f"  {len(all_results)} accessions checked · "
          f"{sum(1 for r in all_results if r[0] == 'ok')} ok · "
          f"{sum(1 for r in all_results if r[0] == 'suspect')} suspect · "
          f"{sum(1 for r in all_results if r[0] == 'error')} error · "
          f"{sum(1 for r in all_results if r[0] == 'unverifiable')} unverifiable")
    if bad:
        print(f"\n  Needs attention:")
        for verdict, panel, acc, declared, reason in bad:
            print(f"    {verdict:<8} {panel:<10} {acc:<18} {declared}")
    print(f"{'=' * 78}\n")

    return 1 if (args.strict and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
