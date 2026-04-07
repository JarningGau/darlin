import csv
import os
import subprocess
import sys
from pathlib import Path


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([os.path.abspath("src"), env.get("PYTHONPATH", "")]).strip(os.pathsep)
    return subprocess.run(
        [sys.executable, "-m", "darlin.cli", *args],
        env=env,
        text=True,
        capture_output=True,
    )


def test_cli_bulk_dna_run_produces_outputs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna/L141_CA_R1.fq.gz")
    fq2 = Path("tests/data/bulkdna/L141_CA_R2.fq.gz")
    assert fq1.exists()
    assert fq2.exists()

    sample_id = "L141_CA"
    outdir = tmp_path / "out"

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        sample_id,
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(outdir),
        "--threads",
        "1",
        "--keep-pear",
        "--sample-n",
        "200",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr

    sample_dir = outdir / sample_id
    assert sample_dir.exists()

    log_file = sample_dir / f"{sample_id}.log"
    assembled_fastq = sample_dir / "pear" / "pear.assembled.fastq"
    extracted = sample_dir / "extracted.tsv"
    filtered = sample_dir / "filtered.tsv"
    combo_dir = sample_dir / "reads_1_u_1_l_0.01"
    denoised_barcodes = combo_dir / "denoised_barcodes.tsv"
    annotated = combo_dir / "annotated.tsv"
    alleles_csv = combo_dir / f"{sample_id}_alleles.csv"

    for p in [log_file, assembled_fastq, extracted, filtered, denoised_barcodes, annotated, alleles_csv]:
        assert p.exists(), f"Expected output missing: {p}"

    with alleles_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames is not None
        assert "md5" in reader.fieldnames
        assert "UMIs" in reader.fieldnames
