from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ComboPaths:
    dir: Path
    denoised_agg_tsv: Path
    denoised_barcodes_tsv: Path
    annotated_tsv: Path
    alleles_tsv: Path

    def ensure_dir(self) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class BulkPaths:
    sample_dir: Path
    pear_dir: Path
    log_file: Path
    assembled_fastq: Path
    extracted_tsv: Path
    filtered_tsv: Path
    pear_log: Path
    pear_out_prefix: Path

    def ensure_dirs(self) -> None:
        self.sample_dir.mkdir(parents=True, exist_ok=True)
        self.pear_dir.mkdir(parents=True, exist_ok=True)

    def combo_dir(self, reads_cutoff: int, umi_ld: int, lb_hd_relative: float) -> Path:
        return self.sample_dir / f"reads_{reads_cutoff}_u_{umi_ld}_l_{lb_hd_relative:.4g}"

    def combo_paths(self, reads_cutoff: int, umi_ld: int, lb_hd_relative: float) -> ComboPaths:
        d = self.combo_dir(reads_cutoff, umi_ld, lb_hd_relative)
        return ComboPaths(
            dir=d,
            denoised_agg_tsv=d / "denoised_agg.tsv",
            denoised_barcodes_tsv=d / "denoised_barcodes.tsv",
            annotated_tsv=d / "annotated.tsv",
            alleles_tsv=d / "alleles_by_umis.tsv",
        )

    def combo_alleles_tsv(self, reads_cutoff: int, umi_ld: int, lb_hd_relative: float) -> Path:
        return self.combo_paths(reads_cutoff, umi_ld, lb_hd_relative).alleles_tsv


def get_bulk_paths(output_dir: str | Path, sample_id: str) -> BulkPaths:
    if not sample_id or "/" in sample_id or "\\" in sample_id:
        raise ValueError(f"sample_id must be a simple name without path separators: {sample_id!r}")
    base = Path(output_dir)
    sample_dir = base / sample_id
    pear_dir = sample_dir / "pear"
    return BulkPaths(
        sample_dir=sample_dir,
        pear_dir=pear_dir,
        log_file=sample_dir / "run.log",
        assembled_fastq=pear_dir / "pear.assembled.fastq",
        extracted_tsv=sample_dir / "extracted.tsv",
        filtered_tsv=sample_dir / "filtered.tsv",
        pear_log=pear_dir / "pear.log",
        pear_out_prefix=pear_dir / "pear",
    )
