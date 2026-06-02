from __future__ import annotations

import logging
import time
from pathlib import Path
from darlin.bulk.logging import setup_logging
from darlin.bulk.paths import get_bulk_paths
from darlin.bulk.steps import (
    step_annotate_and_finalize,
    step_cleanup_pear,
    step_extract,
    step_extract_paired,
    step_filter,
    step_pear,
)


class _StepTimer:
    """Accumulates per-step wall-clock times for a final summary."""

    def __init__(self, logger: logging.Logger) -> None:
        self._logger = logger
        self._records: list[tuple[str, float]] = []
        self._current: tuple[str, float] | None = None
        self._step_idx = 0
        self.total_steps = 0

    def start(self, name: str) -> None:
        self._finish_current()
        self._step_idx += 1
        label = f"[{self._step_idx}/{self.total_steps}]" if self.total_steps else f"[{self._step_idx}]"
        self._logger.info("%s %s", label, name)
        self._current = (name, time.perf_counter())

    def _finish_current(self) -> None:
        if self._current is not None:
            name, t0 = self._current
            self._records.append((name, time.perf_counter() - t0))
            self._current = None

    def summary(self) -> list[tuple[str, float]]:
        self._finish_current()
        return list(self._records)


def resolve_bulk_primers(*, locus: str) -> tuple[int, str, str]:
    # Lazy imports: do not break `--help` when optional deps missing.
    from Bio.Seq import Seq  # type: ignore
    from darlin.bulk.validate import require_valid_locus
    from darlin_core.config.amplicon_configs import load_carlin_config_by_locus  # type: ignore

    require_valid_locus(locus)
    config = load_carlin_config_by_locus(locus=locus)
    unedited_bc_len = len(config.carlin_sequence)
    p5_seq = config.sequence.primer5
    p3_seq = config.sequence.secondary_sequence + config.sequence.primer3
    p3_rc = str(Seq(p3_seq).reverse_complement())
    p5_rc = str(Seq(p5_seq).reverse_complement())
    return unedited_bc_len, p3_rc, p5_rc


def resolve_bulk_primers_paired(*, locus: str) -> tuple[int, str, str]:
    """Forward P3 and P5 strings from darlin_core config for PE85+350 R2 primer matching."""
    from darlin.bulk.validate import require_valid_locus
    from darlin_core.config.amplicon_configs import load_carlin_config_by_locus  # type: ignore

    require_valid_locus(locus)
    config = load_carlin_config_by_locus(locus=locus)
    unedited_bc_len = len(config.carlin_sequence)
    p5_seq = config.sequence.primer5
    p3_seq = config.sequence.secondary_sequence + config.sequence.primer3
    return unedited_bc_len, p3_seq, p5_seq


def _fmt_size(path: Path) -> str:
    """Human-readable file size, or '?' if the file is missing."""
    try:
        n = path.stat().st_size
    except OSError:
        return "?"
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _build_replay_cmd(
    *,
    sample_id: str,
    fq1: str | None,
    fq2: str | None,
    assembled_fq: str | None,
    output_dir: str,
    locus: str,
    protocol: str,
    umi_len: int,
    min_bc_len: int,
    reads_cutoff_list: list[int],
    pear_path: str,
    threads: int,
    log_level: str,
    skip_pear: bool,
    keep_pear: bool,
    test: bool,
    sample_n: int | None,
    umi_ld_list: list[int],
    lb_hd_relative_list: list[float],
    show_progress: bool,
) -> str:
    parts: list[str] = [
        "darlin", "bulk", "run",
        "--sample-id", sample_id,
        "--output-dir", str(output_dir),
        "--locus", locus,
        "--protocol", protocol,
        "--umi-len", str(umi_len),
        "--min-bc-len", str(min_bc_len),
        "--pear-path", pear_path,
        "--threads", str(threads),
        "--log-level", log_level,
    ]
    parts.append("--reads-cutoff")
    parts.extend(str(v) for v in reads_cutoff_list)
    if fq1 is not None:
        parts.extend(["--fq1", fq1])
    if fq2 is not None:
        parts.extend(["--fq2", fq2])
    if skip_pear:
        parts.append("--skip-pear")
        if assembled_fq is not None:
            parts.extend(["--assembled-fq", assembled_fq])
    if keep_pear:
        parts.append("--keep-pear")
    if test:
        parts.append("--test")
    if sample_n is not None:
        parts.extend(["--sample-n", str(sample_n)])
    for v in umi_ld_list:
        parts.extend(["--umi-ld", str(v)])
    for v in lb_hd_relative_list:
        parts.extend(["--lb-hd-relative", str(v)])
    if not show_progress:
        parts.append("--no-progress")
    return " ".join(parts)


def _count_steps(run_pear: bool, n_combos: int) -> int:
    """Return the total number of pipeline steps for progress labels."""
    # PEAR + Extract + Filter + (Denoise + Annotate) * n_combos
    return (1 if run_pear else 0) + 2 + 2 * n_combos


def _log_timing_summary(logger: logging.Logger, records: list[tuple[str, float]]) -> None:
    total = sum(t for _, t in records)
    logger.info("--- Timing summary ---")
    for name, secs in records:
        pct = secs / total * 100 if total > 0 else 0
        logger.info("  %-20s %7.1fs  (%5.1f%%)", name, secs, pct)
    logger.info("  %-20s %7.1fs", "Total", total)


def _log_result_summary(
    logger: logging.Logger,
    combo_allele_paths: list[Path],
    elapsed: float,
) -> None:
    import pandas as pd  # type: ignore

    logger.info("================================")
    logger.info("Pipeline completed in %.1fs", elapsed)
    for p in combo_allele_paths:
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            n_alleles = len(df)
            total_umis = int(df["UMIs"].sum()) if "UMIs" in df.columns else 0
            logger.info("  Unique alleles: %s | Total UMIs: %s | %s", n_alleles, total_umis, p)
    logger.info("================================")


def run_bulk_pipeline(
    *,
    sample_id: str,
    fq1: str | None,
    fq2: str | None,
    assembled_fq: str | None = None,
    output_dir: str,
    locus: str = "Col1a1",
    umi_len: int = 12,
    pear_path: str = "pear",
    threads: int = 8,
    min_bc_len: int = 20,
    reads_cutoff_list: list[int] | None = None,
    denoise_iter: int = 1,
    umi_ld_list: list[int] | None = None,
    lb_hd_relative_list: list[float] | None = None,
    skip_pear: bool = False,
    keep_pear: bool = False,
    log_level: str = "INFO",
    test: bool = False,
    sample_n: int | None = None,
    show_progress: bool = True,
    protocol: str = "pe250",
) -> int:
    if protocol == "pe85-r350" and skip_pear:
        raise ValueError("protocol pe85-r350 cannot be combined with --skip-pear")
    if protocol == "pe250" and not skip_pear and assembled_fq is not None:
        raise ValueError("assembled_fq is only valid when skip_pear=True (protocol pe250)")

    paths = get_bulk_paths(output_dir=output_dir, sample_id=sample_id)
    paths.ensure_dirs()

    level = getattr(logging, log_level.upper(), logging.INFO)
    logger = setup_logging(paths.log_file, level)
    t0 = time.perf_counter()

    if not umi_ld_list:
        umi_ld_list = [1]
    if not lb_hd_relative_list:
        lb_hd_relative_list = [0.01]
    if not reads_cutoff_list:
        reads_cutoff_list = [1]

    max_reads_effective = sample_n if sample_n is not None else (2500 if test else None)

    run_pear = protocol == "pe250" and not skip_pear

    n_combos = len(reads_cutoff_list) * len(umi_ld_list) * len(lb_hd_relative_list)
    timer = _StepTimer(logger)
    timer.total_steps = _count_steps(run_pear=run_pear, n_combos=n_combos)

    replay_cmd = _build_replay_cmd(
        sample_id=sample_id,
        fq1=fq1,
        fq2=fq2,
        assembled_fq=assembled_fq,
        output_dir=output_dir,
        locus=locus,
        protocol=protocol,
        umi_len=umi_len,
        min_bc_len=min_bc_len,
        reads_cutoff_list=reads_cutoff_list,
        pear_path=pear_path,
        threads=threads,
        log_level=log_level,
        skip_pear=skip_pear,
        keep_pear=keep_pear,
        test=test,
        sample_n=sample_n,
        umi_ld_list=umi_ld_list,
        lb_hd_relative_list=lb_hd_relative_list,
        show_progress=show_progress,
    )

    logger.info("--------------------------------")
    logger.info("Starting bulk pipeline for sample: %s", sample_id)
    logger.info("  Protocol: %s", protocol)
    if fq1 is not None and fq2 is not None:
        logger.info("  Input R1: %s (%s)", Path(fq1).name, _fmt_size(Path(fq1)))
        logger.info("  Input R2: %s (%s)", Path(fq2).name, _fmt_size(Path(fq2)))
    elif assembled_fq is not None:
        logger.info("  Assembled: %s (%s)", Path(assembled_fq).name, _fmt_size(Path(assembled_fq)))
    logger.info("  Output:   %s", paths.sample_dir)
    logger.info("  Locus:    %s", locus)
    logger.info("  Threads:  %s", threads)
    if max_reads_effective is not None:
        logger.info("  Max reads: %s (capped)", max_reads_effective)
    logger.debug("Run command (replay): %s", replay_cmd)
    logger.info("--------------------------------")

    if protocol == "pe85-r350":
        unedited_bc_len, p3_fwd, p5_fwd = resolve_bulk_primers_paired(locus=locus)
    else:
        unedited_bc_len, p3_rc, p5_rc = resolve_bulk_primers(locus=locus)

    assembled: Path | None = None
    if run_pear:
        if fq1 is None or fq2 is None:
            raise ValueError("fq1 and fq2 are required when PEAR runs (protocol pe250 without --skip-pear)")
        timer.start("PEAR")
        assembled = step_pear(
            fq1=fq1,
            fq2=fq2,
            paths=paths,
            pear_path=pear_path,
            threads=threads,
            logger=logger,
        )
    elif protocol == "pe250":
        assembled = Path(assembled_fq) if assembled_fq is not None else paths.assembled_fastq
        logger.info("PEAR (skipped)")
        if not assembled.exists():
            raise FileNotFoundError(f"Assembled FASTQ file not found: {assembled}")
    else:
        logger.info("PEAR (skipped): protocol pe85-r350 uses paired FASTQs without PEAR assembly")

    max_reads = max_reads_effective
    timer.start("Extract")
    if protocol == "pe85-r350":
        if fq1 is None or fq2 is None:
            raise ValueError("fq1 and fq2 are required for protocol pe85-r350")
        extracted_tsv = step_extract_paired(
            fq1=fq1,
            fq2=fq2,
            umi_len=umi_len,
            p3_seq=p3_fwd,
            p5_seq=p5_fwd,
            paths=paths,
            max_reads=max_reads,
            logger=logger,
            show_progress=show_progress,
        )
    else:
        if assembled is None:
            raise RuntimeError("internal error: assembled FASTQ path missing for PE250 extract")
        extracted_tsv = step_extract(
            assembled_fastq=assembled,
            umi_len=umi_len,
            p3_seq=p3_rc,
            p5_seq=p5_rc,
            paths=paths,
            max_reads=max_reads,
            logger=logger,
            show_progress=show_progress,
        )

    timer.start("Filter")
    filtered_tsv = step_filter(
        extracted_tsv=extracted_tsv,
        min_bc_len=min_bc_len,
        paths=paths,
        logger=logger,
    )

    combo_allele_paths: list[Path] = []
    for reads_cutoff in reads_cutoff_list:
        for umi_ld in umi_ld_list:
            for lb_rel in lb_hd_relative_list:
                combo = paths.combo_paths(reads_cutoff=reads_cutoff, umi_ld=umi_ld, lb_hd_relative=lb_rel)
                combo.ensure_dir()

                from darlin.bulk.plots import write_post_denoise_plots, write_post_finalize_plots, write_pre_denoise_plots
                from darlin.bulk.steps import step_denoise  # lazy import

                combo_label = f"reads_cutoff={reads_cutoff}, umi_ld={umi_ld}, lb_hd_relative={lb_rel}"
                timer.start(f"Denoise ({combo_label})")
                write_pre_denoise_plots(
                    filtered_tsv=filtered_tsv,
                    out_dir=combo.dir,
                    reads_cutoff=reads_cutoff,
                    unedited_bc_len=unedited_bc_len,
                    logger=logger,
                )
                _denoised_agg_tsv, denoised_barcodes_tsv = step_denoise(
                    filtered_tsv=filtered_tsv,
                    reads_cutoff=reads_cutoff,
                    denoise_iter=denoise_iter,
                    umi_ld=umi_ld,
                    lb_hd_relative=lb_rel,
                    combo=combo,
                    logger=logger,
                    show_progress=show_progress,
                )
                write_post_denoise_plots(
                    denoised_agg_tsv=combo.denoised_agg_tsv,
                    out_dir=combo.dir,
                    unedited_bc_len=unedited_bc_len,
                    logger=logger,
                )

                timer.start("Annotate")
                alleles_tsv = step_annotate_and_finalize(
                    denoised_barcodes_tsv=denoised_barcodes_tsv,
                    locus=locus,
                    min_bc_len=min_bc_len,
                    sample_id=sample_id,
                    combo=combo,
                    logger=logger,
                )
                write_post_finalize_plots(
                    alleles_tsv=alleles_tsv,
                    out_dir=combo.dir,
                    unedited_bc_len=unedited_bc_len,
                    logger=logger,
                )
                combo_allele_paths.append(alleles_tsv)

    step_cleanup_pear(paths=paths, keep_pear=keep_pear, logger=logger)
    elapsed = time.perf_counter() - t0

    _log_timing_summary(logger, timer.summary())
    _log_result_summary(logger, combo_allele_paths, elapsed)
    return 0
