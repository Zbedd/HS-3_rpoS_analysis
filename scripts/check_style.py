"""Audit every drawing module against the house figure style.

    .venv/bin/python scripts/check_style.py           # report
    .venv/bin/python scripts/check_style.py --strict  # non-zero exit on any finding

Six checks, one per way a module can drift out from under `viz`:

    raw-hex         a colour literal outside viz/style.py
    figsize         a figure size that is not a journal column width
    rcparams        rcParams / plt.style set outside viz/style.py
    no-viz          a drawing module that never imports viz
    panel-contract  a panels module that builds a figure or writes a file
    paths           sys.path surgery, or a module counting `..` from itself

A hex literal is the check that matters most: `viz.style` can only be the
single source of truth for colour if no other module names one.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from paths import REPO_ROOT  # noqa: E402

STYLE = Path("src/viz/style.py")

# Diagnostics are exempt from the geometry and colour rules: they are QC
# overlays sized to the plate they draw, read on screen at whatever zoom the
# reader needs, and never typeset. Anything that reaches a manuscript is not
# on this list.
EXEMPT = {
    Path("src/assays/cr/annotate.py"): "per-plate QC overlay, not typeset",
}

# A drawing module is detected by what it does, not by its name: anything that
# touches an axes or pyplot is expected to take its look from viz. `report.py`
# and `promoter.py` name figures in prose and draw none, which a filename list
# would flag and this does not.
DRAWS = re.compile(r"^\s*(ax|fig)[0-9_]*\.|plt\.", re.M)

HEX = re.compile(r"#[0-9a-fA-F]{6}\b")
FIGSIZE = re.compile(r"figsize\s*=\s*(.+)")
RCPARAMS = re.compile(r"rcParams|plt\.style\.")
BAD_PATH = re.compile(r"sys\.path|parent\.parent")
COLUMN = re.compile(r"fs\.COL_|config\.FIG_SIZE|COL_SINGLE|COL_ONEHALF|COL_DOUBLE")

# viz owns geometry and this auditor quotes the patterns it looks for; neither
# is a module drifting away from the style.
SKIP = (Path("src/viz"), Path("scripts/check_style.py"))


def _sources() -> list[Path]:
    out = []
    for base in ("src", "scripts"):
        out += sorted((REPO_ROOT / base).rglob("*.py"))
    rels = [p.relative_to(REPO_ROOT) for p in out]
    return [r for r in rels
            if not any(r == s or s in r.parents for s in SKIP)]


def audit() -> dict[str, list[tuple[Path, int, str]]]:
    found: dict[str, list[tuple[Path, int, str]]] = {
        k: [] for k in ("raw-hex", "figsize", "rcparams", "no-viz",
                        "panel-contract", "paths")}
    for rel in _sources():
        if rel == STYLE:
            continue
        text = (REPO_ROOT / rel).read_text(encoding="utf-8")
        lines = text.splitlines()
        exempt = rel in EXEMPT
        imports_viz = bool(re.search(r"^\s*(import viz|from viz)", text, re.M))
        drawing = bool(DRAWS.search(text)) and rel.parts[0] == "src"

        if drawing and not imports_viz and not exempt:
            found["no-viz"].append((rel, 0, "draws but never imports viz"))

        for n, line in enumerate(lines, 1):
            code = line.split("#")[0] if not HEX.search(line) else line
            if not exempt:
                for hx in HEX.findall(line):
                    found["raw-hex"].append((rel, n, hx))
            m = FIGSIZE.search(code)
            if m and not exempt and not COLUMN.search(m.group(1)):
                found["figsize"].append((rel, n, m.group(1).strip().rstrip(",)")))
            if RCPARAMS.search(code) and "RC_PARAMS" not in code:
                found["rcparams"].append((rel, n, code.strip()[:60]))
            if BAD_PATH.search(code) and rel not in (Path("src/paths.py"),
                                                     Path("scripts/check_style.py")):
                found["paths"].append((rel, n, code.strip()[:60]))
            if rel.name == "panels.py" and re.search(
                    r"plt\.(figure|subplots)\(|\.savefig\(|fs\.save\(", code):
                found["panel-contract"].append((rel, n, code.strip()[:60]))
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if anything is reported")
    args = ap.parse_args()

    found = audit()
    total = sum(len(v) for v in found.values())

    for check, hits in found.items():
        if not hits:
            print(f"[ok]   {check}")
            continue
        by_file: dict[Path, list[tuple[int, str]]] = {}
        for rel, n, what in hits:
            by_file.setdefault(rel, []).append((n, what))
        print(f"[{len(hits):>3}] {check}")
        for rel, items in sorted(by_file.items()):
            shown = ", ".join(f"{w}" for _, w in items[:6])
            more = f" (+{len(items) - 6} more)" if len(items) > 6 else ""
            first = items[0][0]
            where = f"{rel}:{first}" if first else str(rel)
            print(f"         {where}  {shown}{more}")

    if EXEMPT:
        print("\nexempt:")
        for rel, why in EXEMPT.items():
            print(f"  {rel} — {why}")

    print(f"\n{total} finding(s) across {len(_sources())} module(s)")
    return 1 if (args.strict and total) else 0


if __name__ == "__main__":
    raise SystemExit(main())
