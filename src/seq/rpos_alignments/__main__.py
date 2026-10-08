"""Run the manuscript RpoS alignment workflow."""
from __future__ import annotations

import argparse
from pathlib import Path

from .paths import DEFAULT_WORKDIR
from .rpos import run


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--exploratory", action="store_true",
                        help="also render RpoS domain cartoons and sequence logos")
    args = parser.parse_args(argv)
    run(args.workdir, exploratory=args.exploratory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
