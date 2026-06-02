from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from darlin.paths_common import validate_sample_id


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
    final_tsv: Path

    def ensure_dirs(self) -> None:
        self.sample_dir.mkdir(parents=True, exist_ok=True)
        self.diagnostics_dir.mkdir(parents=True, exist_ok=True)


def get_scrna_paths(output_dir: str | Path, sample_id: str) -> ScrnaPaths:
    validate_sample_id(sample_id)
    sample_dir = Path(output_dir) / sample_id
    diagnostics_dir = sample_dir / "diagnostic_plots"
    return ScrnaPaths(
        sample_dir=sample_dir,
        diagnostics_dir=diagnostics_dir,
        log_file=sample_dir / "run.log",
        extracted_tsv=sample_dir / "step1_extracted.tsv",
        denoised_tsv=sample_dir / "step2_denoised.tsv",
        qc_tsv=sample_dir / "step3_qc.tsv",
        cell_summary_tsv=sample_dir / "step3_qc_capture_oligo_carryover_data.tsv",
        annotated_tsv=sample_dir / "step4_annotated.tsv",
        final_tsv=sample_dir / "step4_final.tsv",
    )
