from __future__ import annotations

import argparse
import sys


def add_scmulti_command(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "scmulti",
        help="Recover lineage information from single-cell multi-modal data",
    )
    p.set_defaults(func=_scmulti_main)


def _scmulti_main(_: argparse.Namespace) -> int:
    print(
        "darlin scmulti is not implemented in this repo yet.\n"
        "This command currently exists for backward-compatible CLI shape only.",
        file=sys.stderr,
    )
    return 2
