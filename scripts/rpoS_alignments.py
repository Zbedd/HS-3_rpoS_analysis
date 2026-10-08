#!/usr/bin/env python3
"""Render the four RpoS alignment supplements and their diagnostic tables.

    python scripts/rpoS_alignments.py
    python scripts/rpoS_alignments.py --exploratory  # domain cartoons and logos

Uses the frozen manuscript sequence caches.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from seq.rpos_alignments.paths import DEFAULT_WORKDIR
from seq.rpos_alignments.rpos import run


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--exploratory", action="store_true",
                        help="also render RpoS domain cartoons and sequence logos")
    args = parser.parse_args()
    run(args.workdir, exploratory=args.exploratory)


if __name__ == "__main__":
    main()
