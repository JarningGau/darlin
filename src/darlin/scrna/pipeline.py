from __future__ import annotations

import logging
import time
from pathlib import Path

from darlin.scrna.logging import setup_logging
from darlin.scrna.paths import get_scrna_paths
from darlin.scrna.protocols import ScrnaProtocol
from darlin.scrna.steps import step_annotate, step_denoise, step_extract, step_qc


def resolve_scrna_primers(*, locus: str) -> tuple[int, str, str]:
    from darlin_core.config.amplicon_configs import load_carlin_config_by_locus  # type: ignore

    config = load_carlin_config_by_locus(locus=locus)
    unedited_bc_len = len(config.carlin_sequence)
    p5_seq = config.sequence.primer5
    p3_seq = config.sequence.secondary_sequence + config.sequence.primer3
    return unedited_bc_len, p3_seq, p5_seq


def resolve_max_reads(*, test: bool, sample_n: int | None) -> int | None:
    if sample_n is not None:
        return sample_n
    if test:
        return 2500
    return None


def run_scrna_pipeline(
    *,
    sample_id: str,
    fq1: str,
    fq2: str,
    output_dir: str,
    locus: str,
    protocol: ScrnaProtocol,
    whitelist_path: Path,
    log_level: str,
    test: bool,
    sample_n: int | None,
    show_progress: bool,
    min_bc_len: int,
    umi_ld: int,
    lb_error_rate: float,
    lb_min_hd: int,
    major_fraction_threshold_molecule: float,
    reads_umis_ratio_cutoff: float,
    reads_cutoff: int,
) -> int:
    paths = get_scrna_paths(output_dir=output_dir, sample_id=sample_id)
    paths.ensure_dirs()

    level = getattr(logging, log_level.upper(), logging.INFO)
    logger = setup_logging(paths.log_file, level)
    started = time.perf_counter()

    logger.info('-' * 40)
    logger.info("Starting scrna pipeline for sample: %s", sample_id)
    logger.info("  Input R1: %s", Path(fq1).name)
    logger.info("  Input R2: %s", Path(fq2).name)
    logger.info("  Output:   %s", paths.sample_dir)
    logger.info("  Locus:    %s", locus)
    logger.info("  Protocol: %s", protocol.name)
    logger.info('-' * 40)

    unedited_bc_len, p3_seq, p5_seq = resolve_scrna_primers(locus=locus)
    max_reads = resolve_max_reads(test=test, sample_n=sample_n)

    step_extract(
        fq1=Path(fq1),
        fq2=Path(fq2),
        protocol=protocol,
        p3_seq=p3_seq,
        p5_seq=p5_seq,
        unedited_bc_len=unedited_bc_len,
        paths=paths,
        logger=logger,
        max_reads=max_reads,
        show_progress=show_progress,
    )
    step_denoise(
        extracted_tsv=paths.extracted_tsv,
        whitelist_path=whitelist_path,
        min_bc_len=min_bc_len,
        umi_ld=umi_ld,
        lb_error_rate=lb_error_rate,
        lb_min_hd=lb_min_hd,
        paths=paths,
        logger=logger,
    )
    step_qc(
        denoised_tsv=paths.denoised_tsv,
        major_fraction_threshold_molecule=major_fraction_threshold_molecule,
        reads_umis_ratio_cutoff=reads_umis_ratio_cutoff,
        reads_cutoff=reads_cutoff,
        paths=paths,
        logger=logger,
    )
    step_annotate(
        qc_tsv=paths.qc_tsv,
        locus=locus,
        min_bc_len=min_bc_len,
        paths=paths,
        logger=logger,
    )

    elapsed = time.perf_counter() - started
    logger.info("Pipeline completed in %.2fs", elapsed)
    return 0
