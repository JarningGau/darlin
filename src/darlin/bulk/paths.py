from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BulkPaths:
    sample_dir: Path
    pear_dir: Path
    log_file: Path
    assembled_fastq: Path
    extracted_tsv: Path
    filtered_tsv: Path

    def combo_dir(self, reads_cutoff: int, umi_ld: int, lb_hd_relative: float) -> Path:
        return self.sample_dir / f"reads_{reads_cutoff}_u_{umi_ld}_l_{lb_hd_relative}"

    def combo_alleles_csv(self, reads_cutoff: int, umi_ld: int, lb_hd_relative: float, sample_id: str) -> Path:
        return self.combo_dir(reads_cutoff, umi_ld, lb_hd_relative) / f"{sample_id}_alleles.csv"


def get_bulk_paths(output_dir: str | Path, sample_id: str) -> BulkPaths:
    base = Path(output_dir)
    sample_dir = base / sample_id
    pear_dir = sample_dir / "pear"
    return BulkPaths(
        sample_dir=sample_dir,
        pear_dir=pear_dir,
        log_file=sample_dir / f"{sample_id}.log",
        assembled_fastq=pear_dir / "pear.assembled.fastq",
        extracted_tsv=sample_dir / "extracted.tsv",
        filtered_tsv=sample_dir / "filtered.tsv",
    )

