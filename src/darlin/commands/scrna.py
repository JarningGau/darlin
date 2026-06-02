from __future__ import annotations

import argparse
import sys
from pathlib import Path

from darlin.cli_types import positive_float, positive_int
from darlin.scrna.validate import ScrnaInputError, require_protocol, require_readable_files, resolve_whitelist_path

HELP_FORMATTER = argparse.ArgumentDefaultsHelpFormatter


def add_scrna_command(subparsers: argparse._SubParsersAction) -> None:
    scrna = subparsers.add_parser(
        "scrna",
        help="Recover lineage information from single-cell RNA-seq data",
        formatter_class=HELP_FORMATTER,
    )
    scrna.set_defaults(func=_scrna_main)

    steps = scrna.add_subparsers(dest="scrna_step", metavar="step")
    steps.required = True

    _add_scrna_run(steps)
    _add_scrna_extract(steps)
    _add_scrna_denoise(steps)
    _add_scrna_qc(steps)
    _add_scrna_annotate(steps)


def _scrna_main(_: argparse.Namespace) -> int:
    return 0


def _add_common_scrna_args(
    p: argparse.ArgumentParser,
    *,
    include_fqs: bool,
) -> None:
    p.add_argument("--sample-id", type=str, required=True, help="Sample ID for output directory naming")
    p.add_argument("--output-dir", type=str, default="./output", help="Base output directory")
    p.add_argument("--locus", type=str, default="Col1a1", help="Locus name (darlin-core config key)")
    p.add_argument("--protocol", type=str, default="10xv3", help="Single-cell protocol preset")
    p.add_argument("--whitelist", type=str, default=None, help="Override protocol default whitelist path")
    p.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )
    p.add_argument("--no-progress", action="store_true", help="Disable tqdm progress bars")
    p.add_argument("--test", action="store_true", help="Test mode: only process first ~2500 read pairs")
    p.add_argument("--sample-n", type=positive_int, default=None, help="Process first N read pairs")
    p.add_argument("--min-bc-len", type=positive_int, default=20, help="Minimum lineage barcode length")
    if include_fqs:
        p.add_argument(
            "--fq1",
            type=str,
            required=True,
            help="Path to sequencer Read 1 (R1) FASTQ (protocol determines whether this is the barcode read)",
        )
        p.add_argument(
            "--fq2",
            type=str,
            required=True,
            help="Path to sequencer Read 2 (R2) FASTQ (protocol determines whether this is the barcode read)",
        )


def _add_scrna_run(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("run", help="Run the full scrna pipeline", formatter_class=HELP_FORMATTER)
    _add_common_scrna_args(p, include_fqs=True)
    p.add_argument("--umi-ld", type=positive_int, default=1, help="UMI clustering threshold")
    p.add_argument("--lb-error-rate", type=positive_float, default=0.01, help="Relative lineage barcode error rate")
    p.add_argument("--lb-min-hd", type=positive_int, default=1, help="Minimum lineage barcode Hamming-distance threshold")
    p.add_argument(
        "--major-fraction-threshold-molecule",
        type=positive_float,
        default=0.8,
        help="Minimum reads fraction for the major LR per (CR, UR)",
    )
    p.add_argument(
        "--k-cutoff",
        type=positive_float,
        default=1.0,
        help="Minimum k = reads/UMIs ratio per corrected cell barcode",
    )
    p.add_argument(
        "--reads-cutoff-per-cell",
        type=positive_int,
        default=1,
        help="Minimum total reads per corrected cell barcode after QC",
    )
    p.add_argument(
        "--reads-cutoff-per-molecule",
        type=positive_int,
        default=1,
        help="Minimum reads per molecule retained during denoise",
    )
    p.set_defaults(func=_scrna_run)


def _add_scrna_extract(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser(
        "extract",
        help="Extract lineage barcode, cell barcode, and UMI from FASTQ pairs",
        formatter_class=HELP_FORMATTER,
    )
    _add_common_scrna_args(p, include_fqs=True)
    p.set_defaults(func=_scrna_extract)


def _add_scrna_denoise(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser(
        "denoise",
        help="Denoise cell barcodes, UMIs, and lineage barcodes",
        formatter_class=HELP_FORMATTER,
    )
    _add_common_scrna_args(p, include_fqs=False)
    p.add_argument(
        "--extracted",
        type=str,
        default=None,
        help="Path to step1_extracted.tsv (defaults to <output-dir>/<sample-id>/step1_extracted.tsv)",
    )
    p.add_argument(
        "--reads-cutoff-per-molecule",
        type=positive_int,
        default=1,
        help="Minimum reads per molecule retained during denoise",
    )
    p.add_argument("--umi-ld", type=positive_int, default=1, help="UMI clustering threshold")
    p.add_argument("--lb-error-rate", type=positive_float, default=0.01, help="Relative lineage barcode error rate")
    p.add_argument("--lb-min-hd", type=positive_int, default=1, help="Minimum lineage barcode Hamming-distance threshold")
    p.set_defaults(func=_scrna_denoise)


def _add_scrna_qc(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("qc", help="Run molecular and cellular QC", formatter_class=HELP_FORMATTER)
    _add_common_scrna_args(p, include_fqs=False)
    p.add_argument(
        "--denoised",
        type=str,
        default=None,
        help="Path to step2_denoised.tsv (defaults to <output-dir>/<sample-id>/step2_denoised.tsv)",
    )
    p.add_argument(
        "--major-fraction-threshold-molecule",
        type=positive_float,
        default=0.8,
        help="Minimum reads fraction for the major LR per (CR, UR)",
    )
    p.add_argument(
        "--k-cutoff",
        type=positive_float,
        default=1.0,
        help="Minimum k = reads/UMIs ratio per corrected cell barcode",
    )
    p.add_argument(
        "--reads-cutoff-per-cell",
        type=positive_int,
        default=1,
        help="Minimum total reads per corrected cell barcode after QC",
    )
    p.set_defaults(func=_scrna_qc)


def _add_scrna_annotate(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("annotate", help="Annotate corrected lineage barcodes", formatter_class=HELP_FORMATTER)
    _add_common_scrna_args(p, include_fqs=False)
    p.add_argument(
        "--qc",
        type=str,
        default=None,
        help="Path to step3_qc.tsv (defaults to <output-dir>/<sample-id>/step3_qc.tsv)",
    )
    p.set_defaults(func=_scrna_annotate)


def _resolve_scrna_paths(output_dir: str, sample_id: str):
    from darlin.scrna.paths import get_scrna_paths

    try:
        return get_scrna_paths(output_dir, sample_id)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return None


def _validate_scrna_common(args: argparse.Namespace):
    try:
        protocol = require_protocol(args.protocol)
        whitelist_path = resolve_whitelist_path(protocol, args.whitelist)
    except ScrnaInputError as exc:
        print(exc, file=sys.stderr)
        return None

    paths = _resolve_scrna_paths(args.output_dir, args.sample_id)
    if paths is None:
        return None
    return protocol, whitelist_path, paths


def _validate_scrna_fqs(args: argparse.Namespace) -> bool:
    try:
        require_readable_files(
            [
                ("Forward reads (--fq1)", Path(args.fq1)),
                ("Reverse reads (--fq2)", Path(args.fq2)),
            ]
        )
    except ScrnaInputError as exc:
        print(exc, file=sys.stderr)
        return False
    return True


def _scrna_run(args: argparse.Namespace) -> int:
    if not _validate_scrna_fqs(args):
        return 1
    resolved = _validate_scrna_common(args)
    if resolved is None:
        return 1
    protocol, whitelist_path, _paths = resolved

    from darlin.scrna.pipeline import run_scrna_pipeline

    return run_scrna_pipeline(
        sample_id=args.sample_id,
        fq1=args.fq1,
        fq2=args.fq2,
        output_dir=args.output_dir,
        locus=args.locus,
        protocol=protocol,
        whitelist_path=whitelist_path,
        log_level=args.log_level,
        test=bool(args.test),
        sample_n=args.sample_n,
        show_progress=not bool(args.no_progress),
        min_bc_len=args.min_bc_len,
        umi_ld=args.umi_ld,
        lb_error_rate=args.lb_error_rate,
        lb_min_hd=args.lb_min_hd,
        major_fraction_threshold_molecule=args.major_fraction_threshold_molecule,
        k_cutoff=args.k_cutoff,
        reads_cutoff_per_cell=args.reads_cutoff_per_cell,
        reads_cutoff_per_molecule=args.reads_cutoff_per_molecule,
    )


def _scrna_extract(args: argparse.Namespace) -> int:
    if not _validate_scrna_fqs(args):
        return 1
    resolved = _validate_scrna_common(args)
    if resolved is None:
        return 1
    protocol, _whitelist_path, paths = resolved
    paths.ensure_dirs()

    import logging

    from darlin.scrna.logging import setup_logging
    from darlin.scrna.pipeline import resolve_max_reads, resolve_scrna_primers
    from darlin.scrna.steps import step_extract

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    unedited_bc_len, p3_seq, p5_seq = resolve_scrna_primers(locus=args.locus)
    step_extract(
        fq1=Path(args.fq1),
        fq2=Path(args.fq2),
        protocol=protocol,
        p3_seq=p3_seq,
        p5_seq=p5_seq,
        unedited_bc_len=unedited_bc_len,
        paths=paths,
        logger=logger,
        max_reads=resolve_max_reads(test=bool(args.test), sample_n=args.sample_n),
        show_progress=not bool(args.no_progress),
    )
    return 0


def _scrna_denoise(args: argparse.Namespace) -> int:
    resolved = _validate_scrna_common(args)
    if resolved is None:
        return 1
    _protocol, whitelist_path, paths = resolved

    import logging

    from darlin.scrna.logging import setup_logging
    from darlin.scrna.steps import step_denoise

    extracted_tsv = Path(args.extracted) if args.extracted else paths.extracted_tsv
    try:
        require_readable_files([("Extracted table (--extracted or default output)", extracted_tsv)])
    except ScrnaInputError as exc:
        print(exc, file=sys.stderr)
        print("Run `darlin scrna extract` first, or pass `--extracted`.", file=sys.stderr)
        return 1

    paths.ensure_dirs()
    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_denoise(
        extracted_tsv=extracted_tsv,
        whitelist_path=whitelist_path,
        min_bc_len=args.min_bc_len,
        reads_cutoff_per_molecule=args.reads_cutoff_per_molecule,
        umi_ld=args.umi_ld,
        lb_error_rate=args.lb_error_rate,
        lb_min_hd=args.lb_min_hd,
        paths=paths,
        logger=logger,
    )
    return 0


def _scrna_qc(args: argparse.Namespace) -> int:
    resolved = _validate_scrna_common(args)
    if resolved is None:
        return 1
    _protocol, _whitelist_path, paths = resolved

    import logging

    from darlin.scrna.logging import setup_logging
    from darlin.scrna.steps import step_qc

    denoised_tsv = Path(args.denoised) if args.denoised else paths.denoised_tsv
    try:
        require_readable_files([("Denoised table (--denoised or default output)", denoised_tsv)])
    except ScrnaInputError as exc:
        print(exc, file=sys.stderr)
        print("Run `darlin scrna denoise` first, or pass `--denoised`.", file=sys.stderr)
        return 1

    paths.ensure_dirs()
    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_qc(
        denoised_tsv=denoised_tsv,
        major_fraction_threshold_molecule=args.major_fraction_threshold_molecule,
        k_cutoff=args.k_cutoff,
        reads_cutoff_per_cell=args.reads_cutoff_per_cell,
        paths=paths,
        logger=logger,
    )
    return 0


def _scrna_annotate(args: argparse.Namespace) -> int:
    resolved = _validate_scrna_common(args)
    if resolved is None:
        return 1
    _protocol, _whitelist_path, paths = resolved

    import logging

    from darlin.scrna.logging import setup_logging
    from darlin.scrna.steps import step_annotate

    qc_tsv = Path(args.qc) if args.qc else paths.qc_tsv
    try:
        require_readable_files([("QC table (--qc or default output)", qc_tsv)])
    except ScrnaInputError as exc:
        print(exc, file=sys.stderr)
        print("Run `darlin scrna qc` first, or pass `--qc`.", file=sys.stderr)
        return 1

    paths.ensure_dirs()
    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_annotate(
        qc_tsv=qc_tsv,
        locus=args.locus,
        min_bc_len=args.min_bc_len,
        paths=paths,
        logger=logger,
    )
    return 0
