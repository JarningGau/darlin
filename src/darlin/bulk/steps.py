from __future__ import annotations

import hashlib
import logging
import shutil
from pathlib import Path

from darlin.bulk.denoise import correct_lineage_and_umi
from darlin.bulk.io import iter_fastq_raw, open_fastq_text
from darlin.bulk.matching import find_all_matches, get_mm_dist
from darlin.bulk.paths import BulkPaths
from darlin.bulk.pear import assemble_pe_reads


def step_pear(
    *,
    fq1: str,
    fq2: str,
    paths: BulkPaths,
    pear_path: str = "pear",
    threads: int = 8,
    logger: logging.Logger,
) -> Path:
    paths.pear_dir.mkdir(parents=True, exist_ok=True)
    pear_out_prefix = paths.pear_dir / "pear"
    pear_log = paths.pear_dir / "pear.log"
    assemble_pe_reads(
        fq1,
        fq2,
        out_prefix=pear_out_prefix,
        log_file=pear_log,
        pear_path=pear_path,
        threads=threads,
        logger=logger,
    )
    if not paths.assembled_fastq.exists():
        raise FileNotFoundError(f"PEAR assembled FASTQ not found: {paths.assembled_fastq}")
    return paths.assembled_fastq


def step_extract(
    *,
    assembled_fastq: str | Path,
    umi_len: int,
    p3_seq: str,
    p5_seq: str,
    paths: BulkPaths,
    max_reads: int | None,
    logger: logging.Logger,
) -> Path:
    import pandas as pd  # type: ignore
    from tqdm import tqdm  # type: ignore

    assembled_fastq = Path(assembled_fastq)
    if not assembled_fastq.exists():
        raise FileNotFoundError(f"Assembled FASTQ not found: {assembled_fastq}")

    p3_mm = get_mm_dist(p3_seq)
    p5_mm = get_mm_dist(p5_seq)

    results: list[tuple[str, str]] = []
    read_count = 0
    with open_fastq_text(str(assembled_fastq)) as fq_handle:
        for (_read_id, seq, _qual) in tqdm(
            iter_fastq_raw(fq_handle),
            desc="Processing reads",
            unit_scale=True,
            unit=" reads",
        ):
            if max_reads is not None and read_count >= max_reads:
                break
            read_count += 1
            umi = seq[:umi_len]
            if "N" in umi:
                continue
            match_result = find_all_matches(seq, p3_seq, p5_seq, p3_mm, p5_mm)
            if len(match_result.p3_matches) == 1 and len(match_result.p5_matches) == 1:
                s = match_result.p3_matches[0].end
                e = match_result.p5_matches[0].start
                lineage_bc = seq[s:e]
                results.append((lineage_bc, umi))

    df = pd.DataFrame(results, columns=["lineage_bc", "UMI"])
    df["bc_len"] = df["lineage_bc"].str.len()

    paths.sample_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = paths.extracted_tsv
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Extracted table saved to: {out_tsv}")
    return out_tsv


def step_filter(
    *,
    extracted_tsv: str | Path,
    min_bc_len: int,
    paths: BulkPaths,
    logger: logging.Logger,
) -> Path:
    import pandas as pd  # type: ignore

    extracted_tsv = Path(extracted_tsv)
    if not extracted_tsv.exists():
        raise FileNotFoundError(f"Extracted table not found: {extracted_tsv}")

    df = pd.read_csv(extracted_tsv, sep="\t")
    df["bc_len"] = df["lineage_bc"].astype(str).str.len()
    df = df[df["bc_len"] >= int(min_bc_len)]

    df = df.groupby(["lineage_bc", "UMI"]).size().reset_index(name="reads")
    df["bc_len"] = df["lineage_bc"].astype(str).str.len()
    df.sort_values(by="reads", ascending=False, inplace=True)

    out_tsv = paths.filtered_tsv
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Filtered table saved to: {out_tsv}")
    return out_tsv


def step_denoise(
    *,
    filtered_tsv: str | Path,
    reads_cutoff: int,
    denoise_iter: int,
    umi_ld: int,
    lb_hd_relative: float,
    paths: BulkPaths,
    logger: logging.Logger,
) -> tuple[Path, Path]:
    import pandas as pd  # type: ignore

    filtered_tsv = Path(filtered_tsv)
    if not filtered_tsv.exists():
        raise FileNotFoundError(f"Filtered table not found: {filtered_tsv}")

    df = pd.read_csv(filtered_tsv, sep="\t")
    if "reads" not in df.columns:
        raise ValueError(f"Expected a 'reads' column in {filtered_tsv}")
    n_before = len(df)
    df = df[df["reads"] >= int(reads_cutoff)].copy()
    if df.empty:
        raise ValueError(
            f"No rows remain after applying reads_cutoff={reads_cutoff} "
            f"(from {n_before} aggregated pairs in {filtered_tsv})."
        )
    # correct_lineage_and_umi expects (lineage_bc, UMI, reads)
    agg, _mapping, stats = correct_lineage_and_umi(
        df.rename(columns={"reads": "n_reads"}),
        umi_col="UMI",
        bc_col="lineage_bc",
        count_col="n_reads",
        n_iter=int(denoise_iter),
        umi_ld=int(umi_ld),
        lb_hd_relative=float(lb_hd_relative),
        logger=logger,
    )
    logger.info(f"Denoising stats (umi_ld={umi_ld}, lb_hd_relative={lb_hd_relative}): {stats}")

    combo_dir = paths.combo_dir(reads_cutoff=int(reads_cutoff), umi_ld=umi_ld, lb_hd_relative=lb_hd_relative)
    combo_dir.mkdir(parents=True, exist_ok=True)

    out_agg = combo_dir / "denoised_agg.tsv"
    agg.to_csv(out_agg, sep="\t", index=False)

    from Bio.Seq import Seq  # type: ignore

    agg2 = agg.groupby("lineage_bc_corr").size().reset_index(name="UMIs")
    agg2.sort_values(by="UMIs", ascending=False, inplace=True)
    agg2["query"] = [str(Seq(s).reverse_complement()) for s in agg2["lineage_bc_corr"].astype(str)]
    out_bc = combo_dir / "denoised_barcodes.tsv"
    agg2.to_csv(out_bc, sep="\t", index=False)

    logger.info(f"Denoised outputs saved to: {out_agg} and {out_bc}")
    return out_agg, out_bc


def step_annotate(
    *,
    denoised_barcodes_tsv: str | Path,
    locus: str,
    min_bc_len: int,
    paths: BulkPaths,
    reads_cutoff: int,
    umi_ld: int,
    lb_hd_relative: float,
    logger: logging.Logger,
) -> Path:
    import pandas as pd  # type: ignore
    from Bio.Seq import Seq  # type: ignore
    from darlinpy import analyze_sequences  # type: ignore

    denoised_barcodes_tsv = Path(denoised_barcodes_tsv)
    if not denoised_barcodes_tsv.exists():
        raise FileNotFoundError(f"Denoised barcode table not found: {denoised_barcodes_tsv}")

    agg2 = pd.read_csv(denoised_barcodes_tsv, sep="\t")
    if "query" not in agg2.columns:
        sequences = agg2["lineage_bc_corr"].astype(str).tolist()
        sequences_rc = [str(Seq(s).reverse_complement()) for s in sequences]
        agg2["query"] = sequences_rc
        agg2.to_csv(denoised_barcodes_tsv, sep="\t", index=False)
        logger.info(f"Added `query` column and updated: {denoised_barcodes_tsv}")
    else:
        sequences_rc = agg2["query"].astype(str).tolist()

    results_allele = analyze_sequences(
        sequences_rc,
        config=locus,
        min_sequence_length=int(min_bc_len),
        verbose=False,
    ).to_df()

    combo_dir = paths.combo_dir(reads_cutoff=reads_cutoff, umi_ld=umi_ld, lb_hd_relative=lb_hd_relative)
    combo_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = combo_dir / "annotated.tsv"
    results_allele.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Annotation results saved to: {out_tsv}")
    return out_tsv


def _concat_and_md5(aligned_query: str, aligned_ref: str) -> str:
    concat_str = str(aligned_query) + str(aligned_ref)
    return hashlib.md5(concat_str.encode("utf-8")).hexdigest()


def step_finalize(
    *,
    denoised_barcodes_tsv: str | Path,
    annotated_tsv: str | Path,
    sample_id: str,
    paths: BulkPaths,
    reads_cutoff: int,
    umi_ld: int,
    lb_hd_relative: float,
    logger: logging.Logger,
) -> Path:
    import pandas as pd  # type: ignore

    denoised_barcodes_tsv = Path(denoised_barcodes_tsv)
    annotated_tsv = Path(annotated_tsv)
    if not denoised_barcodes_tsv.exists():
        raise FileNotFoundError(f"Denoised barcode table not found: {denoised_barcodes_tsv}")
    if not annotated_tsv.exists():
        raise FileNotFoundError(f"Annotated table not found: {annotated_tsv}")

    agg2 = pd.read_csv(denoised_barcodes_tsv, sep="\t")
    results_allele = pd.read_csv(annotated_tsv, sep="\t")

    if "query" not in agg2.columns:
        raise ValueError(
            "Expected `query` column in denoised barcodes table. "
            "Use `denoised_barcodes.tsv` from `darlin bulk denoise` (or after `annotate` has added `query`)."
        )

    final = agg2.merge(results_allele, on="query", how="left")
    keep_cols = ["query", "UMIs", "mutations", "confidence", "aligned_query", "aligned_ref"]
    final = final[[c for c in keep_cols if c in final.columns]].copy()

    final["md5"] = final.apply(lambda r: _concat_and_md5(r["aligned_query"], r["aligned_ref"]), axis=1)

    final2 = (
        final.drop(columns=["query"], errors="ignore")
        .groupby("md5", as_index=False)
        .agg(
            {
                "UMIs": "sum",
                "mutations": "first",
                "confidence": "first",
                "aligned_query": "first",
                "aligned_ref": "first",
            }
        )
        .sort_values(by="UMIs", ascending=False)
        .reset_index(drop=True)
    )

    combo_dir = paths.combo_dir(reads_cutoff=reads_cutoff, umi_ld=umi_ld, lb_hd_relative=lb_hd_relative)
    combo_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = paths.combo_alleles_tsv(reads_cutoff=reads_cutoff, umi_ld=umi_ld, lb_hd_relative=lb_hd_relative)
    final2.to_csv(out_tsv, sep="\t", index=False)
    logger.info(f"Final alleles saved to: {out_tsv}")
    return out_tsv


def step_cleanup_pear(*, paths: BulkPaths, keep_pear: bool, logger: logging.Logger) -> None:
    if keep_pear:
        return
    if paths.pear_dir.exists():
        logger.info("Cleaning up PEAR output directory...")
        shutil.rmtree(paths.pear_dir)
        logger.info(f"PEAR output directory removed: {paths.pear_dir}")

