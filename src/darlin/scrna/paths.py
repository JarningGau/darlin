from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ScrnaPaths:
    sample_dir: Path
    diagnostics_dir: Path
    log_file: Path
    extracted_tsv: Path
    denoised_tsv: Path
    qc_tsv: Path
    cell_summary_tsv: Path
    annotated_tsv: Path

    def ensure_dirs(self) -> None:
        self.sample_dir.mkdir(parents=True, exist_ok=True)
        self.diagnostics_dir.mkdir(parents=True, exist_ok=True)


def get_scrna_paths(output_dir: str | Path, sample_id: str) -> ScrnaPaths:
    if "/" in sample_id or "\\" in sample_id:
        raise ValueError("sample_id must be a single path segment without path separators")
    sample_dir = Path(output_dir) / sample_id
    diagnostics_dir = sample_dir / "diagnostics"
    return ScrnaPaths(
        sample_dir=sample_dir,
        diagnostics_dir=diagnostics_dir,
        log_file=sample_dir / "run.log",
        extracted_tsv=sample_dir / "extracted.tsv",
        denoised_tsv=sample_dir / "denoised.tsv",
        qc_tsv=sample_dir / "qc.tsv",
        cell_summary_tsv=sample_dir / "cell_summary.tsv",
        annotated_tsv=sample_dir / "annotated.tsv",
    )

