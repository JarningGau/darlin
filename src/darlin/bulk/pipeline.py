from __future__ import annotations

import logging
from pathlib import Path

from darlin.bulk.logging import setup_logging
from darlin.bulk.paths import get_bulk_paths
from darlin.bulk.steps import (
    step_annotate,
    step_cleanup_pear,
    step_extract,
    step_filter,
    step_finalize,
    step_pear,
)


def resolve_bulk_primers(*, locus: str) -> tuple[int, str, str]:
    # Lazy imports: do not break `--help` when optional deps missing.
    from Bio.Seq import Seq  # type: ignore
    from darlinpy.config.amplicon_configs import load_carlin_config_by_locus  # type: ignore

    config = load_carlin_config_by_locus(locus=locus)
    unedited_bc_len = len(config.carlin_sequence)
    p5_seq = config.sequence.primer5
    p3_seq = config.sequence.secondary_sequence + config.sequence.primer3
    p3_rc = str(Seq(p3_seq).reverse_complement())
    p5_rc = str(Seq(p5_seq).reverse_complement())
    return unedited_bc_len, p3_rc, p5_rc


def run_bulk_pipeline(
    *,
    sample_id: str,
    fq1: str,
    fq2: str,
    output_dir: str,
    locus: str = "Col1a1",
    umi_len: int = 12,
    pear_path: str = "pear",
    threads: int = 8,
    min_bc_len: int = 20,
    reads_cutoff: int = 1,
    denoise_iter: int = 1,
    umi_ld_list: list[int] | None = None,
    lb_hd_relative_list: list[float] | None = None,
    skip_pear: bool = False,
    keep_pear: bool = False,
    log_level: str = "INFO",
    test: bool = False,
    sample_n: int | None = None,
) -> int:
    paths = get_bulk_paths(output_dir=output_dir, sample_id=sample_id)
    paths.sample_dir.mkdir(parents=True, exist_ok=True)
    paths.pear_dir.mkdir(parents=True, exist_ok=True)

    level = getattr(logging, log_level.upper(), logging.INFO)
    logger = setup_logging(paths.log_file, level)

    logger.info("--------------------------------")
    logger.info(f"Starting bulk pipeline for sample: {sample_id}")
    logger.info(f"Output directory: {paths.sample_dir}")
    logger.info("--------------------------------")

    if not umi_ld_list:
        umi_ld_list = [1]
    if not lb_hd_relative_list:
        lb_hd_relative_list = [0.01]

    _unedited_bc_len, p3_rc, p5_rc = resolve_bulk_primers(locus=locus)

    if not skip_pear:
        logger.info("#### Step: PEAR")
        assembled = step_pear(
            fq1=fq1,
            fq2=fq2,
            paths=paths,
            pear_path=pear_path,
            threads=threads,
            logger=logger,
        )
    else:
        assembled = paths.assembled_fastq
        logger.info("#### Step: PEAR (skipped)")
        if not assembled.exists():
            raise FileNotFoundError(f"Assembled FASTQ file not found: {assembled}")

    max_reads = sample_n if sample_n is not None else (2500 if test else None)
    logger.info("#### Step: Extract")
    extracted_tsv = step_extract(
        assembled_fastq=assembled,
        umi_len=umi_len,
        p3_seq=p3_rc,
        p5_seq=p5_rc,
        paths=paths,
        max_reads=max_reads,
        logger=logger,
    )

    logger.info("#### Step: Filter")
    filtered_tsv = step_filter(
        extracted_tsv=extracted_tsv,
        min_bc_len=min_bc_len,
        reads_cutoff=reads_cutoff,
        paths=paths,
        logger=logger,
    )

    # For now, we keep the original behavior: loop over parameter combinations and
    # write final alleles per-combo.
    for umi_ld in umi_ld_list:
        for lb_rel in lb_hd_relative_list:
            combo_dir = paths.combo_dir(reads_cutoff=reads_cutoff, umi_ld=umi_ld, lb_hd_relative=lb_rel)
            combo_dir.mkdir(parents=True, exist_ok=True)

            # Denoise step currently writes into combo_dir; annotate/finalize consume those outputs.
            from darlin.bulk.steps import step_denoise  # lazy import

            logger.info(f"#### Combo: reads_cutoff={reads_cutoff}, umi_ld={umi_ld}, lb_hd_relative={lb_rel}")
            _denoised_agg_tsv, denoised_barcodes_tsv = step_denoise(
                filtered_tsv=filtered_tsv,
                reads_cutoff=reads_cutoff,
                denoise_iter=denoise_iter,
                umi_ld=umi_ld,
                lb_hd_relative=lb_rel,
                paths=paths,
                logger=logger,
            )

            logger.info("#### Step: Annotate")
            annotated_tsv = step_annotate(
                denoised_barcodes_tsv=denoised_barcodes_tsv,
                locus=locus,
                min_bc_len=min_bc_len,
                paths=paths,
                reads_cutoff=reads_cutoff,
                umi_ld=umi_ld,
                lb_hd_relative=lb_rel,
                logger=logger,
            )

            logger.info("#### Step: Finalize")
            _alleles_csv = step_finalize(
                denoised_barcodes_tsv=denoised_barcodes_tsv,
                annotated_tsv=annotated_tsv,
                sample_id=sample_id,
                paths=paths,
                reads_cutoff=reads_cutoff,
                umi_ld=umi_ld,
                lb_hd_relative=lb_rel,
                logger=logger,
            )

    step_cleanup_pear(paths=paths, keep_pear=keep_pear, logger=logger)
    logger.info("Processing completed successfully!")
    return 0

