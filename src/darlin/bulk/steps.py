from __future__ import annotations

import hashlib
import logging
import shutil
from pathlib import Path

from darlin.bulk.denoise import correct_lineage_and_umi
from darlin.bulk.io import iter_fastq_paired, iter_fastq_raw, open_fastq_text
from darlin.bulk.matching import find_all_matches, get_mm_dist
from darlin.bulk.paths import BulkPaths, ComboPaths
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
    assemble_pe_reads(
        fq1,
        fq2,
        out_prefix=paths.pear_out_prefix,
        log_file=paths.pear_log,
        pear_path=pear_path,
        threads=threads,
        logger=logger,
        log_path_bases=(Path.cwd(), paths.sample_dir),
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
    show_progress: bool = True,
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
    skipped_umi_n = 0
    skipped_no_dual_match = 0
    with open_fastq_text(str(assembled_fastq)) as fq_handle:
        for (_read_id, seq, _qual) in tqdm(
            iter_fastq_raw(fq_handle),
            desc="Processing reads",
            unit_scale=True,
            unit=" reads",
            disable=not show_progress,
        ):
            if max_reads is not None and read_count >= max_reads:
                break
            read_count += 1
            umi = seq[:umi_len]
            if "N" in umi:
                skipped_umi_n += 1
                continue
            match_result = find_all_matches(seq, p3_seq, p5_seq, p3_mm, p5_mm)
            if len(match_result.p3_matches) == 1 and len(match_result.p5_matches) == 1:
                s = match_result.p3_matches[0].end
                e = match_result.p5_matches[0].start
                lineage_bc = seq[s:e]
                results.append((lineage_bc, umi))
            else:
                skipped_no_dual_match += 1

    df = pd.DataFrame(results, columns=["LB", "UB"])
    df["LB_len"] = df["LB"].str.len()

    out_tsv = paths.extracted_tsv
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(
        "Extract: reads_scanned=%s rows_written=%s skipped_umi_with_N=%s skipped_no_unique_primer_pair=%s -> %s",
        read_count,
        len(results),
        skipped_umi_n,
        skipped_no_dual_match,
        out_tsv,
    )
    return out_tsv


def step_extract_paired(
    *,
    fq1: str | Path,
    fq2: str | Path,
    umi_len: int,
    p3_seq: str,
    p5_seq: str,
    paths: BulkPaths,
    max_reads: int | None,
    logger: logging.Logger,
    show_progress: bool = True,
) -> Path:
    """
    PE85+350 paired extraction: UMI from R1 5'; lineage barcode on R2 between P5 and P3 (forward primers).
    """
    import pandas as pd  # type: ignore
    from tqdm import tqdm  # type: ignore

    fq1 = Path(fq1)
    fq2 = Path(fq2)
    if not fq1.exists():
        raise FileNotFoundError(f"R1 FASTQ not found: {fq1}")
    if not fq2.exists():
        raise FileNotFoundError(f"R2 FASTQ not found: {fq2}")

    p3_mm = get_mm_dist(p3_seq)
    p5_mm = get_mm_dist(p5_seq)

    results: list[tuple[str, str]] = []
    read_count = 0
    skipped_umi_n = 0
    skipped_no_dual_match = 0

    with open_fastq_text(str(fq1)) as h1, open_fastq_text(str(fq2)) as h2:
        for (_id1, seq1, _q1, _id2, seq2, _q2) in tqdm(
            iter_fastq_paired(h1, h2),
            desc="Processing reads",
            unit_scale=True,
            unit=" reads",
            disable=not show_progress,
        ):
            if max_reads is not None and read_count >= max_reads:
                break
            read_count += 1
            umi = seq1[:umi_len]
            if "N" in umi:
                skipped_umi_n += 1
                continue
            match_result = find_all_matches(seq2, p3_seq, p5_seq, p3_mm, p5_mm)
            if len(match_result.p3_matches) == 1 and len(match_result.p5_matches) == 1:
                s = match_result.p5_matches[0].end
                e = match_result.p3_matches[0].start
                lineage_bc = seq2[s:e]
                results.append((lineage_bc, umi))
            else:
                skipped_no_dual_match += 1

    df = pd.DataFrame(results, columns=["LB", "UB"])
    df["LB_len"] = df["LB"].str.len()

    out_tsv = paths.extracted_tsv
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.info(
        "Extract (paired): reads_scanned=%s rows_written=%s skipped_umi_with_N=%s "
        "skipped_no_unique_primer_pair=%s -> %s",
        read_count,
        len(results),
        skipped_umi_n,
        skipped_no_dual_match,
        out_tsv,
    )
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
    df["LB_len"] = df["LB"].astype(str).str.len()
    df = df[df["LB_len"] >= int(min_bc_len)]

    df = df.groupby(["LB", "UB"]).size().reset_index(name="reads")
    df["LB_len"] = df["LB"].astype(str).str.len()
    df.sort_values(by="reads", ascending=False, inplace=True)

    out_tsv = paths.filtered_tsv
    df.to_csv(out_tsv, sep="\t", index=False)
    n_rows = len(df)
    sum_reads = int(df["reads"].sum()) if n_rows > 0 else 0
    logger.info("Filter: aggregated_rows=%s sum_reads=%s -> %s", n_rows, sum_reads, out_tsv)
    return out_tsv


def step_denoise(
    *,
    filtered_tsv: str | Path,
    reads_cutoff: int,
    denoise_iter: int,
    umi_ld: int,
    lb_hd_relative: float,
    combo: ComboPaths,
    logger: logging.Logger,
    show_progress: bool = True,
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
    # correct_lineage_and_umi expects (LB, UB, reads)
    agg, _mapping, stats = correct_lineage_and_umi(
        df,
        umi_col="UB",
        bc_col="LB",
        count_col="reads",
        n_iter=int(denoise_iter),
        umi_ld=int(umi_ld),
        lb_hd_relative=float(lb_hd_relative),
        logger=logger,
        show_progress=show_progress,
    )
    agg.to_csv(combo.denoised_agg_tsv, sep="\t", index=False)

    from Bio.Seq import Seq  # type: ignore

    agg2 = agg.groupby("LR").size().reset_index(name="UMIs")
    agg2.sort_values(by="UMIs", ascending=False, inplace=True)
    agg2["query"] = [str(Seq(s).reverse_complement()) for s in agg2["LR"].astype(str)]
    agg2.to_csv(combo.denoised_barcodes_tsv, sep="\t", index=False)

    logger.info(
        "Denoise: umi_ld=%s lb_hd_relative=%s n_input_rows=%s "
        "pairs_before=%s pairs_after=%s umi_merges=%s barcode_merges=%s "
        "-> %s",
        umi_ld,
        lb_hd_relative,
        stats["n_input_rows"],
        stats["n_unique_pairs_before"],
        stats["n_unique_pairs_after"],
        stats["umi_merges"],
        stats["barcode_merges"],
        combo.denoised_barcodes_tsv,
    )
    return combo.denoised_agg_tsv, combo.denoised_barcodes_tsv


def step_annotate(
    *,
    denoised_barcodes_tsv: str | Path,
    locus: str,
    min_bc_len: int,
    combo: ComboPaths,
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
        sequences = agg2["LR"].astype(str).tolist()
        sequences_rc = [str(Seq(s).reverse_complement()) for s in sequences]
        agg2["query"] = sequences_rc
        agg2.to_csv(denoised_barcodes_tsv, sep="\t", index=False)
        logger.info(f"Added `query` column and updated: {denoised_barcodes_tsv}")
    else:
        sequences_rc = agg2["query"].astype(str).tolist()

    n_queries = len(sequences_rc)
    results_allele = analyze_sequences(
        sequences_rc,
        config=locus,
        min_sequence_length=int(min_bc_len),
        verbose=False,
    ).to_df()

    results_allele.to_csv(combo.annotated_tsv, sep="\t", index=False)
    logger.info(
        "Annotation: analyzed_queries=%s rows_out=%s -> %s",
        n_queries,
        len(results_allele),
        combo.annotated_tsv,
    )
    return combo.annotated_tsv


def step_annotate_and_finalize(
    *,
    denoised_barcodes_tsv: str | Path,
    locus: str,
    min_bc_len: int,
    sample_id: str,
    combo: ComboPaths,
    logger: logging.Logger,
) -> Path:
    annotated_tsv = step_annotate(
        denoised_barcodes_tsv=denoised_barcodes_tsv,
        locus=locus,
        min_bc_len=min_bc_len,
        combo=combo,
        logger=logger,
    )
    return step_finalize(
        denoised_barcodes_tsv=denoised_barcodes_tsv,
        annotated_tsv=annotated_tsv,
        sample_id=sample_id,
        combo=combo,
        logger=logger,
    )


def _concat_and_md5(aligned_query: str, aligned_ref: str) -> str:
    concat_str = str(aligned_query) + str(aligned_ref)
    return hashlib.md5(concat_str.encode("utf-8")).hexdigest()


def step_finalize(
    *,
    denoised_barcodes_tsv: str | Path,
    annotated_tsv: str | Path,
    sample_id: str,
    combo: ComboPaths,
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
    input_queries = len(final)
    annotated_mask = final["aligned_query"].notna() & final["aligned_ref"].notna()
    dropped_unannotated = int((~annotated_mask).sum())
    final = final[annotated_mask].copy()
    keep_cols = ["query", "UMIs", "mutations", "confidence", "aligned_query", "aligned_ref"]
    final = final[[c for c in keep_cols if c in final.columns]].copy()

    if final.empty:
        final2 = final.drop(columns=["query"], errors="ignore").copy()
        final2["md5"] = pd.Series(dtype=str)
        final2 = final2.reindex(columns=["md5", "UMIs", "mutations", "confidence", "aligned_query", "aligned_ref"])
    else:
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

    final2.to_csv(combo.alleles_tsv, sep="\t", index=False)
    logger.info(
        "Finalize: input_queries=%s annotated_queries=%s unannotated_queries_dropped=%s final_alleles=%s -> %s",
        input_queries,
        input_queries - dropped_unannotated,
        dropped_unannotated,
        len(final2),
        combo.alleles_tsv,
    )
    return combo.alleles_tsv


def step_cleanup_pear(*, paths: BulkPaths, keep_pear: bool, logger: logging.Logger) -> None:
    if keep_pear:
        return
    if paths.pear_dir.exists():
        logger.info("Cleaning up PEAR output directory...")
        shutil.rmtree(paths.pear_dir)
        logger.info(f"PEAR output directory removed: {paths.pear_dir}")
