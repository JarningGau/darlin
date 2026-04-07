from __future__ import annotations

import argparse
import sys
from pathlib import Path

from darlin.bulk.validate import BulkInputError, require_pear_executable, require_readable_files


def add_bulk_command(subparsers: argparse._SubParsersAction) -> None:
    bulk = subparsers.add_parser(
        "bulk",
        help="Recover lineage information from bulk DNA/RNA data",
    )
    bulk.set_defaults(func=_bulk_main)

    steps = bulk.add_subparsers(dest="bulk_step", metavar="step")
    steps.required = True

    # Step commands will be wired to real implementations in `darlin.bulk.*`.
    _add_bulk_run(steps)
    _add_bulk_pear(steps)
    _add_bulk_extract(steps)
    _add_bulk_filter(steps)
    _add_bulk_denoise(steps)
    _add_bulk_annotate(steps)
    _add_bulk_finalize(steps)


def _bulk_main(_: argparse.Namespace) -> int:
    # If user runs `darlin bulk` without a step, argparse will handle help/usage.
    return 0


def _add_bulk_run(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("run", help="Run the full bulk pipeline")
    _add_common_bulk_args(p, include_fqs=True)
    p.add_argument("--skip-pear", action="store_true", help="Skip PEAR assembly (use existing assembled file)")
    p.add_argument("--keep-pear", action="store_true", help="Keep PEAR output directory after completion")
    p.add_argument("--denoise-iter", type=int, default=1, help="Number of denoising iterations")
    p.add_argument("--umi-ld", type=int, nargs="+", default=[1], help="List of UMI clustering thresholds")
    p.add_argument("--lb-hd-relative", type=float, nargs="+", default=[0.01], help="List of relative barcode HD thresholds")
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.add_argument("--test", action="store_true", help="Test mode: only process first ~2500 reads")
    p.add_argument("--sample-n", type=int, default=None, help="Sample first N reads")
    p.set_defaults(func=_bulk_run)


def _add_bulk_pear(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("pear", help="Assemble paired-end reads (PEAR)")
    _add_common_bulk_args(p, include_fqs=True)
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.set_defaults(func=_bulk_pear)


def _add_bulk_extract(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("extract", help="Extract lineage barcode + UMI from assembled FASTQ")
    _add_common_bulk_args(p, include_fqs=False)
    p.add_argument(
        "--assembled-fq",
        type=str,
        default=None,
        help="Path to assembled FASTQ (defaults to <output-dir>/<sample-id>/pear/pear.assembled.fastq)",
    )
    p.add_argument("--skip-pear", action="store_true", help="Skip PEAR assembly (requires assembled FASTQ to exist)")
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.add_argument("--test", action="store_true", help="Test mode: only process first ~2500 reads")
    p.add_argument("--sample-n", type=int, default=None, help="Sample first N reads")
    p.set_defaults(func=_bulk_extract)


def _add_bulk_filter(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("filter", help="Filter/aggregate extracted reads")
    _add_common_bulk_args(p, include_fqs=False, include_reads_cutoff=False)
    p.add_argument(
        "--extracted",
        type=str,
        default=None,
        help="Path to extracted.tsv (defaults to <output-dir>/<sample-id>/extracted.tsv)",
    )
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.set_defaults(func=_bulk_filter)


def _add_bulk_denoise(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("denoise", help="Denoise lineage barcodes and UMIs")
    _add_common_bulk_args(p, include_fqs=False)
    p.add_argument(
        "--filtered",
        type=str,
        default=None,
        help="Path to filtered.tsv (defaults to <output-dir>/<sample-id>/filtered.tsv)",
    )
    p.add_argument("--denoise-iter", type=int, default=1, help="Number of denoising iterations")
    p.add_argument("--umi-ld", type=int, default=1, help="UMI clustering threshold")
    p.add_argument("--lb-hd-relative", type=float, default=0.01, help="Relative barcode HD threshold")
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.set_defaults(func=_bulk_denoise)


def _add_bulk_annotate(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("annotate", help="Annotate alleles using darlinpy")
    _add_common_bulk_args(p, include_fqs=False)
    p.add_argument("--umi-ld", type=int, default=1, help="UMI clustering threshold (for output dir naming)")
    p.add_argument("--lb-hd-relative", type=float, default=0.01, help="Relative barcode HD threshold (for output dir naming)")
    p.add_argument(
        "--denoised-barcodes",
        type=str,
        required=True,
        help="Path to denoised_barcodes.tsv (typically produced by `darlin bulk denoise`)",
    )
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.set_defaults(func=_bulk_annotate)


def _add_bulk_finalize(steps: argparse._SubParsersAction) -> None:
    p = steps.add_parser("finalize", help="Finalize outputs (CSV + diagnostics)")
    _add_common_bulk_args(p, include_fqs=False)
    p.add_argument("--umi-ld", type=int, default=1, help="UMI clustering threshold (for output dir naming)")
    p.add_argument("--lb-hd-relative", type=float, default=0.01, help="Relative barcode HD threshold (for output dir naming)")
    p.add_argument("--denoised-barcodes", type=str, required=True, help="Path to denoised_barcodes.tsv")
    p.add_argument("--annotated", type=str, required=True, help="Path to annotated.tsv")
    p.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level")
    p.set_defaults(func=_bulk_finalize)


def _add_common_bulk_args(
    p: argparse.ArgumentParser,
    *,
    include_fqs: bool,
    include_reads_cutoff: bool = True,
) -> None:
    p.add_argument("--sample-id", type=str, required=True, help="Sample ID for output directory naming")
    p.add_argument("--output-dir", type=str, default="./output", help="Base output directory")
    p.add_argument("--locus", type=str, default="Col1a1", help="Locus name (darlinpy config key)")
    p.add_argument("--umi-len", type=int, default=12, help="UMI length")
    p.add_argument("--min-bc-len", type=int, default=20, help="Minimum barcode length")
    if include_reads_cutoff:
        p.add_argument(
            "--reads-cutoff",
            type=int,
            default=1,
            help="Minimum read support per (lineage barcode, UMI) pair; applied at denoising (not at filter)",
        )
    p.add_argument("--pear-path", type=str, default="pear", help="Path to PEAR executable")
    p.add_argument("--threads", type=int, default=8, help="Number of threads for PEAR")
    if include_fqs:
        p.add_argument("--fq1", type=str, required=True, help="Path to forward reads FASTQ file")
        p.add_argument("--fq2", type=str, required=True, help="Path to reverse reads FASTQ file")


def _bulk_run(args: argparse.Namespace) -> int:
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.pipeline import run_bulk_pipeline

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    try:
        if args.skip_pear:
            require_readable_files(
                [("Assembled FASTQ (default: <output-dir>/<sample-id>/pear/pear.assembled.fastq)", paths.assembled_fastq)]
            )
        else:
            require_readable_files(
                [
                    ("Forward reads (--fq1)", args.fq1),
                    ("Reverse reads (--fq2)", args.fq2),
                ]
            )
            require_pear_executable(args.pear_path)
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    return run_bulk_pipeline(
        sample_id=args.sample_id,
        fq1=args.fq1,
        fq2=args.fq2,
        output_dir=args.output_dir,
        locus=args.locus,
        umi_len=args.umi_len,
        pear_path=args.pear_path,
        threads=args.threads,
        min_bc_len=args.min_bc_len,
        reads_cutoff=args.reads_cutoff,
        denoise_iter=args.denoise_iter,
        umi_ld_list=list(args.umi_ld),
        lb_hd_relative_list=list(args.lb_hd_relative),
        skip_pear=bool(args.skip_pear),
        keep_pear=bool(args.keep_pear),
        log_level=args.log_level,
        test=bool(args.test),
        sample_n=args.sample_n,
    )


def _bulk_pear(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.steps import step_pear

    try:
        require_readable_files(
            [
                ("Forward reads (--fq1)", args.fq1),
                ("Reverse reads (--fq2)", args.fq2),
            ]
        )
        require_pear_executable(args.pear_path)
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_pear(fq1=args.fq1, fq2=args.fq2, paths=paths, pear_path=args.pear_path, threads=args.threads, logger=logger)
    return 0


def _bulk_extract(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.pipeline import resolve_bulk_primers
    from darlin.bulk.steps import step_extract

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    assembled = Path(args.assembled_fq) if args.assembled_fq else paths.assembled_fastq
    try:
        require_readable_files([("Assembled FASTQ (--assembled-fq or default pear output)", assembled)])
    except BulkInputError as e:
        print(e, file=sys.stderr)
        print("Run `darlin bulk pear` first, or pass `--assembled-fq`.", file=sys.stderr)
        return 1

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))

    _unedited, p3_rc, p5_rc = resolve_bulk_primers(locus=args.locus)

    max_reads = args.sample_n if args.sample_n is not None else (2500 if args.test else None)
    step_extract(
        assembled_fastq=assembled,
        umi_len=args.umi_len,
        p3_seq=p3_rc,
        p5_seq=p5_rc,
        paths=paths,
        max_reads=max_reads,
        logger=logger,
    )
    return 0


def _bulk_filter(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.steps import step_filter

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    extracted = Path(args.extracted) if args.extracted else paths.extracted_tsv
    try:
        require_readable_files([("Extracted TSV (--extracted or default)", extracted)])
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_filter(
        extracted_tsv=extracted,
        min_bc_len=args.min_bc_len,
        paths=paths,
        logger=logger,
    )
    return 0


def _bulk_denoise(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.steps import step_denoise

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    filtered = Path(args.filtered) if args.filtered else paths.filtered_tsv
    try:
        require_readable_files([("Filtered TSV (--filtered or default)", filtered)])
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_denoise(
        filtered_tsv=filtered,
        reads_cutoff=args.reads_cutoff,
        denoise_iter=args.denoise_iter,
        umi_ld=args.umi_ld,
        lb_hd_relative=args.lb_hd_relative,
        paths=paths,
        logger=logger,
    )
    return 0


def _bulk_annotate(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.steps import step_annotate

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    try:
        require_readable_files([("Denoised barcodes (--denoised-barcodes)", args.denoised_barcodes)])
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_annotate(
        denoised_barcodes_tsv=args.denoised_barcodes,
        locus=args.locus,
        min_bc_len=args.min_bc_len,
        paths=paths,
        reads_cutoff=args.reads_cutoff,
        umi_ld=args.umi_ld,
        lb_hd_relative=args.lb_hd_relative,
        logger=logger,
    )
    return 0


def _bulk_finalize(args: argparse.Namespace) -> int:
    import logging

    from darlin.bulk.logging import setup_logging
    from darlin.bulk.paths import get_bulk_paths
    from darlin.bulk.steps import step_finalize

    paths = get_bulk_paths(args.output_dir, args.sample_id)
    try:
        require_readable_files(
            [
                ("Denoised barcodes (--denoised-barcodes)", args.denoised_barcodes),
                ("Annotated TSV (--annotated)", args.annotated),
            ]
        )
    except BulkInputError as e:
        print(e, file=sys.stderr)
        return 1

    logger = setup_logging(paths.log_file, getattr(logging, args.log_level.upper(), logging.INFO))
    step_finalize(
        denoised_barcodes_tsv=args.denoised_barcodes,
        annotated_tsv=args.annotated,
        sample_id=args.sample_id,
        paths=paths,
        reads_cutoff=args.reads_cutoff,
        umi_ld=args.umi_ld,
        lb_hd_relative=args.lb_hd_relative,
        logger=logger,
    )
    return 0
