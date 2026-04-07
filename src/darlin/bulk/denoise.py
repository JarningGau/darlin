from __future__ import annotations

import logging
from collections import Counter
from typing import Optional, Sequence, Union


def hamming_dist(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def neighbors_hd1(s: str) -> list[str]:
    bases = ("A", "C", "G", "T")
    out: list[str] = []
    for i, ch in enumerate(s):
        for b in bases:
            if b != ch:
                out.append(s[:i] + b + s[i + 1 :])
    return out


def collapse_within_hd(items, max_hd: int):
    counts = dict(items)
    seqs = sorted(counts, key=lambda s: counts[s], reverse=True)
    parent = {s: s for s in seqs}

    for i, s in enumerate(seqs):
        if parent[s] != s:
            continue
        c_hi = counts[s]
        for t in seqs[i + 1 :]:
            if parent[t] != t:
                continue
            if len(t) != len(s):
                continue
            if hamming_dist(s, t) <= max_hd:
                c_lo = counts[t]
                if c_hi >= 2 * c_lo - 1:
                    parent[t] = s
    return parent


def _to_df(data, bc_col: str, umi_col: str, count_col: Optional[str]):
    import pandas as pd  # type: ignore

    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        df = pd.DataFrame(data)

    if bc_col not in df.columns or umi_col not in df.columns:
        if df.shape[1] >= 2 and bc_col not in df.columns and umi_col not in df.columns:
            cols = list(df.columns)
            rename_map = {}
            if len(cols) >= 1:
                rename_map[cols[0]] = bc_col
            if len(cols) >= 2:
                rename_map[cols[1]] = umi_col
            if len(cols) >= 3 and count_col and count_col not in df.columns:
                rename_map[cols[2]] = count_col
            df = df.rename(columns=rename_map)

    missing = [c for c in [bc_col, umi_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing required column(s): {missing}. Available: {list(df.columns)}")

    for c in [bc_col, umi_col]:
        df[c] = df[c].astype(str).str.upper()
    df = df.dropna(subset=[bc_col, umi_col])

    if count_col and count_col in df.columns:
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(1).astype(int)
    else:
        df["__count__"] = 1
        count_col = "__count__"

    return df, count_col


def correct_lineage_and_umi(
    data: Union["pd.DataFrame", Sequence],  # noqa: F821
    umi_col: str = "UB",
    bc_col: str = "LB",
    count_col: Optional[str] = None,
    n_iter: int = 2,
    umi_ld: int = 2,
    lb_hd_relative: float = 0.01,
    logger=None,
):
    import pandas as pd  # type: ignore
    from tqdm import tqdm  # type: ignore
    from umi_tools import UMIClusterer  # type: ignore

    if logger is None:
        logger = logging.getLogger(__name__)

    df, count_col = _to_df(data, bc_col=bc_col, umi_col=umi_col, count_col=count_col)
    out = df[[bc_col, umi_col, count_col]].copy()
    out["UR"] = out[umi_col]
    out["LR"] = out[bc_col]

    total_umi_merges = 0
    total_bc_merges = 0

    clusterer = UMIClusterer(cluster_method="directional")

    for i in range(n_iter):
        logger.info(f"Iteration {i+1}/{n_iter}")
        out["__bc_len__"] = out["LR"].str.len().astype(int)
        bc_parent_total = {}
        for blen, sub in tqdm(
            out.groupby(["__bc_len__"]),
            desc="Collapsing barcodes (length-aware HD global)",
            leave=True,
        ):
            cnt = Counter(dict(sub.groupby("LR")[count_col].sum()))
            if isinstance(blen, (tuple, list)):
                blen = int(blen[0])
            else:
                blen = int(blen)
            if blen == 0:
                continue
            hd_thresh_len = max(int(round(lb_hd_relative * blen)), 1)
            parent = collapse_within_hd(cnt.items(), max_hd=hd_thresh_len)
            bc_parent_total.update(parent)
        before = out["LR"].ne(out["LR"].map(lambda b: bc_parent_total.get(b, b))).sum()
        out["LR"] = out["LR"].map(lambda b: bc_parent_total.get(b, b))
        total_bc_merges += int(before)

        parent_all = {}
        for bc_val, sub in tqdm(
            out.groupby("LR"),
            desc="Collapsing UMIs with umi_tools",
            leave=True,
        ):
            _ = bc_val
            cnt_series = sub.groupby("UR")[count_col].sum()
            if cnt_series.empty:
                continue
            umi_counts_bytes = {umi.encode(): int(c) for umi, c in cnt_series.items()}
            umi_groups = clusterer(umi_counts_bytes, threshold=umi_ld)

            for group in umi_groups:
                representative_bytes = max(group, key=lambda u: umi_counts_bytes.get(u, 0))
                representative_str = representative_bytes.decode()
                for umi_bytes in group:
                    umi_str = umi_bytes.decode()
                    parent_all[umi_str] = representative_str

        before = out["UR"].ne(out["UR"].map(lambda u: parent_all.get(u, u))).sum()
        out["UR"] = out["UR"].map(lambda u: parent_all.get(u, u))
        total_umi_merges += int(before)

    agg = (
        out.groupby(["LR", "UR"], as_index=False)[count_col]
        .sum()
        .rename(columns={count_col: "reads"})
    )

    mapping = out[[bc_col, umi_col, "LR", "UR"]].copy().drop_duplicates()

    stats = {
        "n_input_rows": int(len(df)),
        "n_unique_pairs_before": int(df.groupby([bc_col, umi_col]).size().shape[0]),
        "n_unique_pairs_after": int(agg.shape[0]),
        "umi_merges": total_umi_merges,
        "barcode_merges": total_bc_merges,
    }
    return agg, mapping, stats

