from __future__ import annotations

import argparse
import sys


def add_scrna_command(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "scrna",
        help="Recover lineage information from single-cell RNA-seq data",
    )
    p.set_defaults(func=_scrna_main)


def _scrna_main(_: argparse.Namespace) -> int:
    print(
        "darlin scrna is not implemented in this repo yet.\n"
        "This command currently exists for backward-compatible CLI shape only.",
        file=sys.stderr,
    )
    return 2
