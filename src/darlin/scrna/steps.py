from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd  # type: ignore
from tqdm import tqdm  # type: ignore

from darlin.scrna.io import iter_fastq_paired, load_whitelist, open_text_maybe_gzip
from darlin.scrna.logging import (
    log_count,
    log_count_with_pct,
    log_molecule_summary,
    log_output,
    log_step_title,
)
from darlin.scrna.matching import find_all_matches, get_mm_dist
from darlin.scrna.paths import ScrnaPaths
from darlin.scrna.plots import write_extract_plots, write_qc_plots
from darlin.scrna.protocols import ScrnaProtocol


def step_extract(
    *,
    fq1: Path,
    fq2: Path,
    protocol: ScrnaProtocol,
    p3_seq: str,
    p5_seq: str,
    unedited_bc_len: int,
    paths: ScrnaPaths,
    logger: logging.Logger,
    max_reads: int | None,
    show_progress: bool,
) -> pd.DataFrame:
    p3_mm = get_mm_dist(p3_seq)
    p5_mm = get_mm_dist(p5_seq)

    rows: list[tuple[str, str, str, int]] = []
    total_reads = 0
    matched_reads = 0
    skipped_barcode_n = 0
    skipped_no_dual_match = 0

    with open_text_maybe_gzip(fq1) as handle1, open_text_maybe_gzip(fq2) as handle2:
        iterator = iter_fastq_paired(handle1, handle2)
        iterator = tqdm(
            iterator,
            desc="Processing reads",
            unit=" reads",
            unit_scale=True,
            disable=not show_progress,
        )
        for _id1, seq1, _qual1, _id2, seq2, _qual2 in iterator:
            total_reads += 1
            if max_reads is not None and total_reads > max_reads:
                break

            reads = {"fq1": seq1, "fq2": seq2}
            barcode_seq = reads[protocol.barcode_read]
            darlin_seq = reads[protocol.darlin_read]

            cb = barcode_seq[: protocol.cb_len]
            ub = barcode_seq[protocol.cb_len : protocol.cb_len + protocol.umi_len]
            if "N" in cb or "N" in ub:
                skipped_barcode_n += 1
                continue

            match_result = find_all_matches(darlin_seq, p3_seq, p5_seq, p3_mm, p5_mm)
            if len(match_result.p3_matches) == 1 and len(match_result.p5_matches) == 1:
                start = match_result.p5_matches[0].end
                end = match_result.p3_matches[0].start
                lb = darlin_seq[start:end]
                rows.append((lb, cb, ub, len(lb)))
                matched_reads += 1
            else:
                skipped_no_dual_match += 1

    df = pd.DataFrame(rows, columns=["LB", "CB", "UB", "LB_len"])
    df.to_csv(paths.extracted_tsv, sep="\t", index=False)
    write_extract_plots(df, paths.diagnostics_dir, unedited_bc_len=unedited_bc_len)
    reads_scanned = total_reads if max_reads is None else min(total_reads, max_reads)
    log_step_title(logger, "Extract")
    log_count(logger, "Total reads", reads_scanned)
    log_count_with_pct(logger, "Matched reads", matched_reads, reads_scanned)
    log_count_with_pct(logger, "Skipped (barcode with N)", skipped_barcode_n, reads_scanned)
    log_count_with_pct(
        logger,
        "Skipped (no matched flanking sequence)",
        skipped_no_dual_match,
        reads_scanned,
    )
    log_output(logger, paths.extracted_tsv)
    return df


def _neighbors_hd1(seq: str) -> list[str]:
    bases = ("A", "C", "G", "T")
    out: list[str] = []
    for i, ch in enumerate(seq):
        for base in bases:
            if base != ch:
                out.append(seq[:i] + base + seq[i + 1 :])
    return out


def _correct_cb_to_whitelist(values: pd.Series, whitelist: set[str]) -> list[str | None]:
    corrected: list[str | None] = []
    for cb in values.astype(str):
        if cb in whitelist:
            corrected.append(cb)
            continue
        hits = [candidate for candidate in _neighbors_hd1(cb) if candidate in whitelist]
        corrected.append(hits[0] if len(hits) == 1 else None)
    return corrected


def _correct_umis_per_cell(
    df: pd.DataFrame,
    *,
    cb_col: str,
    umi_col: str,
    count_col: str,
    threshold: int,
) -> pd.DataFrame:
    from umi_tools import UMIClusterer  # type: ignore

    clusterer = UMIClusterer(cluster_method="directional")
    corrected: list[pd.DataFrame] = []

    for cr, sub in df.groupby(cb_col, sort=False):
        counts = sub.groupby(umi_col)[count_col].sum()
        umi_counts = {umi.encode(): int(n) for umi, n in counts.items()}
        groups = clusterer(umi_counts, threshold=threshold)

        umi_to_rep: dict[str, str] = {}
        for group in groups:
            rep = max(group, key=lambda umi: umi_counts[umi]).decode()
            for umi in group:
                umi_to_rep[umi.decode()] = rep

        sub = sub.copy()
        sub["UR"] = sub[umi_col].map(lambda umi: umi_to_rep.get(str(umi), str(umi)))
        corrected.append(sub)

    if not corrected:
        df = df.copy()
        df["UR"] = pd.Series(dtype=str)
        return df
    return pd.concat(corrected, ignore_index=True)


def _hamming_dist(a: str, b: str) -> int | None:
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))


def _collapse_within_hd(items: list[tuple[str, int]], *, max_hd: int) -> dict[str, str]:
    ordered = sorted(items, key=lambda item: (-item[1], item[0]))
    representatives: list[str] = []
    mapping: dict[str, str] = {}

    for seq, _count in ordered:
        assigned = None
        for rep in representatives:
            hd = _hamming_dist(seq, rep)
            if hd is not None and hd <= max_hd:
                assigned = rep
                break
        if assigned is None:
            representatives.append(seq)
            mapping[seq] = seq
        else:
            mapping[seq] = assigned
    return mapping


def _correct_lb_per_cell(
    df: pd.DataFrame,
    *,
    cr_col: str,
    lb_col: str,
    lb_len_col: str,
    count_col: str,
    error_rate: float,
    min_hd: int,
) -> pd.DataFrame:
    corrected: list[pd.DataFrame] = []

    for (_cr, lb_len), sub in df.groupby([cr_col, lb_len_col], sort=False):
        counts = sub.groupby(lb_col)[count_col].sum()
        if counts.empty:
            corrected.append(sub)
            continue
        hd_threshold = max(int(round(error_rate * int(lb_len))), min_hd)
        lb_to_rep = _collapse_within_hd(list(counts.items()), max_hd=hd_threshold)
        sub = sub.copy()
        sub["LR"] = sub[lb_col].map(lambda lb: lb_to_rep.get(str(lb), str(lb)))
        corrected.append(sub)

    if not corrected:
        df = df.copy()
        df["LR"] = pd.Series(dtype=str)
        return df
    return pd.concat(corrected, ignore_index=True)


def step_denoise(
    *,
    extracted_tsv: Path,
    whitelist_path: Path,
    min_bc_len: int,
    umi_ld: int,
    lb_error_rate: float,
    lb_min_hd: int,
    paths: ScrnaPaths,
    logger: logging.Logger,
) -> pd.DataFrame:
    df = pd.read_csv(extracted_tsv, sep="\t")
    df = (
        df
        .groupby(['LB', 'CB', 'UB', 'LB_len'], as_index=False)
        .size()
        .rename(columns={"size": "reads"})
    )
    df = df[df["LB_len"] >= min_bc_len].copy()
    log_step_title(logger, "Before denoise")
    log_molecule_summary(logger, df)
    log_step_title(logger, "Denoise (correct sequencing errors)")
    ## CB -> CR
    whitelist = load_whitelist(whitelist_path)
    df["CR"] = _correct_cb_to_whitelist(df["CB"], whitelist)
    df = df[df["CR"].notna()].copy()
    df = (
        df
        .groupby(['LB', 'CR', 'UB', 'LB_len'], as_index=False)
        .agg(reads=('reads', 'sum'))
    )
    ## UR -> UR
    df = _correct_umis_per_cell(df, cb_col="CR", umi_col="UB", count_col="reads", threshold=umi_ld)
    ## LR -> LR
    df = _correct_lb_per_cell(
        df,
        cr_col="CR",
        lb_col="LB",
        lb_len_col="LB_len",
        count_col="reads",
        error_rate=lb_error_rate,
        min_hd=lb_min_hd,
    )
    df = (
        df
        .groupby(['LR', 'CR', 'UR', 'LB_len'], as_index=False)
        .agg(reads=('reads', 'sum'))
    )
    df.to_csv(paths.denoised_tsv, sep="\t", index=False)
    log_step_title(logger, "After sequencing error correction")
    log_molecule_summary(logger, df)
    log_output(logger, paths.denoised_tsv)
    return df


def step_qc(
    *,
    denoised_tsv: Path,
    major_fraction_threshold_molecule: float,
    reads_umis_ratio_cutoff: float,
    reads_cutoff: int,
    paths: ScrnaPaths,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    log_step_title(logger, "QC (remove PCR artifacts)")
    df = pd.read_csv(denoised_tsv, sep="\t")
    df["group_reads"] = df.groupby(["CR", "UR"])["reads"].transform("sum")
    df["reads_fraction"] = df["reads"] / df["group_reads"]
    df_major = df[df["reads_fraction"] >= major_fraction_threshold_molecule].copy()

    cell_summary = (
        df_major.groupby("CR")
        .agg(n_reads=("reads", "sum"), n_UR=("UR", "nunique"))
        .reset_index()
    )
    if len(cell_summary) > 0:
        cell_summary["k"] = cell_summary["n_reads"] / cell_summary["n_UR"]
    else:
        cell_summary["k"] = pd.Series(dtype=float)

    df_merged = df_major.merge(cell_summary[["CR", "k"]], on="CR", how="left")
    df_after_k = df_merged[df_merged["k"] >= reads_umis_ratio_cutoff].copy()
    df_final = df_after_k[df_after_k["reads"] >= reads_cutoff].copy()
    if len(df_final) > 0:
        df_final["n_LR"] = df_final.groupby("CR")["LR"].transform("nunique")
    else:
        df_final["n_LR"] = pd.Series(dtype=int)

    df_final.to_csv(paths.qc_tsv, sep="\t", index=False)
    cell_summary.to_csv(paths.cell_summary_tsv, sep="\t", index=False)
    write_qc_plots(
        df,
        cell_summary,
        df_final,
        paths.diagnostics_dir,
        major_fraction_threshold_molecule=major_fraction_threshold_molecule,
    )
    reads_with_amplification_error = int(df["reads"].sum() - df_major["reads"].sum())
    log_step_title(logger, "After PCR chimera removal")
    log_molecule_summary(logger, df_major)
    log_count(logger, "Reads with amplification error", reads_with_amplification_error)

    reads_removed_coca = int(df_major["reads"].sum() - df_final["reads"].sum())
    log_step_title(logger, "After capture oligo carryover removal")
    log_molecule_summary(logger, df_final)
    log_count(logger, "Reads removed as capture-oligo carryover", reads_removed_coca)
    log_output(logger, paths.qc_tsv)
    return df_final, cell_summary


def _validate_grouped_annotation_uniqueness(source: pd.DataFrame) -> pd.DataFrame:
    mapping = source[["LR", "mutation", "aligned_LR", "aligned_ref"]].drop_duplicates()
    counts = mapping.groupby("LR", dropna=False).size()
    bad_lrs = counts[counts > 1].index.tolist()
    if bad_lrs:
        preview = ", ".join(str(lr) for lr in bad_lrs[:5])
        raise ValueError(f"Non-unique annotation mapping for LR: {preview}")
    return mapping


def step_annotate(
    *,
    qc_tsv: Path,
    locus: str,
    min_bc_len: int,
    paths: ScrnaPaths,
    logger: logging.Logger,
) -> pd.DataFrame:
    from darlin_core import analyze_sequences  # type: ignore

    def _write_grouped_counts(final_df: pd.DataFrame) -> pd.DataFrame:
        # grouped_cols = ["n_UMIs", "CR", "LR", "mutation", "aligned_LR", "aligned_ref"]
        grouped_cols = ["n_UMIs", "CR", "mutation", "aligned_LR", "aligned_ref"]
        source = final_df.rename(columns={"mutations": "mutation"})
        # annotation_map = _validate_grouped_annotation_uniqueness(source)
        # [note]
        # LR -> mutation: many to one due to sequence or PCR errors
        # darlin_core.analyze_sequences() is able to correct these errors
        grouped = (
            source.groupby(["CR", "mutation", "aligned_LR", "aligned_ref"], dropna=False, as_index=False)["UR"]
            .nunique()
            .rename(columns={"UR": "n_UMIs"})
        )
        # grouped = grouped.merge(annotation_map, on="LR", how="left")
        grouped = grouped[grouped_cols]
        grouped.to_csv(paths.final_tsv, sep="\t", index=False)
        return grouped

    df = pd.read_csv(qc_tsv, sep="\t")
    query = df["LR"].dropna().astype(str).drop_duplicates().tolist()
    if not query:
        raise ValueError(f"No lineage barcodes found in QC output: {qc_tsv}")

    results = analyze_sequences(
        query,
        config=locus,
        min_sequence_length=min_bc_len,
        verbose=False,
    ).to_df()
    merged = df.merge(results, left_on="LR", right_on="query", how="left")
    final_df = merged[["CR", "LR", "UR", "reads", "mutations", "aligned_query", "aligned_ref"]].copy()
    final_df = final_df.rename(columns={"aligned_query": "aligned_LR"})
    final_df.to_csv(paths.annotated_tsv, sep="\t", index=False)
    grouped = _write_grouped_counts(final_df)
    log_step_title(logger, "Annotate")
    log_output(logger, paths.annotated_tsv)
    log_output(logger, paths.final_tsv)
    return final_df
