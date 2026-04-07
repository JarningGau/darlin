from __future__ import annotations

import argparse
import sys

from darlin.commands.bulk import add_bulk_command
from darlin.commands.scrna import add_scrna_command
from darlin.commands.scmulti import add_scmulti_command


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="darlin",
        description=(
            "DARLIN: unified CLI for bulk DNA/RNA, single-cell RNA-seq, "
            "and single-cell multi-modal datasets."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", metavar="command")
    subparsers.required = True

    add_bulk_command(subparsers)
    add_scrna_command(subparsers)
    add_scmulti_command(subparsers)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        return int(args.func(args))
    except BrokenPipeError:
        # allow piping into `head` etc.
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
