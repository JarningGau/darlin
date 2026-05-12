import csv
import io
import logging
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

from darlin.bulk.paths import ComboPaths
from darlin.bulk.steps import step_finalize


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([os.path.abspath("src"), env.get("PYTHONPATH", "")]).strip(os.pathsep)
    return subprocess.run(
        [sys.executable, "-m", "darlin.cli", *args],
        env=env,
        text=True,
        capture_output=True,
    )


def test_cli_bulk_dna_run_pe85_r350_produces_outputs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna-f85r350/C126_CA_R1.fq.gz")
    fq2 = Path("tests/data/bulkdna-f85r350/C126_CA_R2.fq.gz")
    assert fq1.exists()
    assert fq2.exists()

    sample_id = "C126_CA"
    outdir = tmp_path / "out"

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        sample_id,
        "--protocol",
        "pe85-r350",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(outdir),
        "--threads",
        "1",
        "--sample-n",
        "200",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr

    sample_dir = outdir / sample_id
    assert sample_dir.exists()

    log_file = sample_dir / "run.log"
    extracted = sample_dir / "extracted.tsv"
    filtered = sample_dir / "filtered.tsv"
    combo_dir = sample_dir / "reads_1_u_1_l_0.01"
    denoised_barcodes = combo_dir / "denoised_barcodes.tsv"
    annotated = combo_dir / "annotated.tsv"
    alleles_tsv = combo_dir / "alleles_by_umis.tsv"

    pear_assembled = sample_dir / "pear" / "pear.assembled.fastq"
    assert not pear_assembled.exists()

    for p in [log_file, extracted, filtered, denoised_barcodes, annotated, alleles_tsv]:
        assert p.exists(), f"Expected output missing: {p}"

    db = pd.read_csv(denoised_barcodes, sep="\t")
    assert (db["query"].astype(str) == db["LR"].astype(str)).all()

    log_text = log_file.read_text()
    assert "Protocol: pe85-r350" in log_text
    assert "Extract (paired):" in log_text
    assert "PEAR (skipped): protocol pe85-r350" in log_text


def test_cli_bulk_run_reads_cutoff_list_produces_multiple_combo_dirs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna-f85r350/C126_CA_R1.fq.gz")
    fq2 = Path("tests/data/bulkdna-f85r350/C126_CA_R2.fq.gz")
    assert fq1.exists()
    assert fq2.exists()

    sample_id = "C126_CA"
    outdir = tmp_path / "out"

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        sample_id,
        "--protocol",
        "pe85-r350",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(outdir),
        "--threads",
        "1",
        "--reads-cutoff",
        "1",
        "2",
        "--sample-n",
        "200",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr

    sample_dir = outdir / sample_id
    assert (sample_dir / "reads_1_u_1_l_0.01" / "alleles_by_umis.tsv").exists()
    assert (sample_dir / "reads_2_u_1_l_0.01" / "alleles_by_umis.tsv").exists()


def test_cli_bulk_run_pe85_r350_rejects_skip_pear(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna-f85r350/C126_CA_R1.fq.gz")
    fq2 = Path("tests/data/bulkdna-f85r350/C126_CA_R2.fq.gz")
    assert fq1.exists()
    assert fq2.exists()

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        "C126_CA",
        "--protocol",
        "pe85-r350",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--skip-pear",
        "--threads",
        "1",
    )
    assert r.returncode == 1
    assert "Traceback" not in r.stderr


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

    log_file = sample_dir / "run.log"
    assembled_fastq = sample_dir / "pear" / "pear.assembled.fastq"
    extracted = sample_dir / "extracted.tsv"
    filtered = sample_dir / "filtered.tsv"
    combo_dir = sample_dir / "reads_1_u_1_l_0.01"
    denoised_barcodes = combo_dir / "denoised_barcodes.tsv"
    annotated = combo_dir / "annotated.tsv"
    alleles_tsv = combo_dir / "alleles_by_umis.tsv"

    for p in [log_file, assembled_fastq, extracted, filtered, denoised_barcodes, annotated, alleles_tsv]:
        assert p.exists(), f"Expected output missing: {p}"

    db = pd.read_csv(denoised_barcodes, sep="\t")
    assert (db["query"].astype(str) == db["LR"].astype(str)).all()

    with alleles_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames is not None
        assert list(reader.fieldnames) == ["LR", "UMIs", "mutation", "aligned_query", "aligned_ref"]

    log_text = log_file.read_text()
    assert "Starting bulk pipeline for sample: L141_CA" in log_text
    assert "Extract: reads_scanned=" in log_text and "rows_written=" in log_text
    assert "Filter: aggregated_rows=" in log_text and "sum_reads=" in log_text
    assert "] Denoise" in log_text
    assert "Denoise: umi_ld=" in log_text and "barcode_merges=" in log_text
    assert "Annotation: analyzed_queries=" in log_text
    assert "Timing summary" in log_text
    assert "Pipeline completed in" in log_text


def test_cli_bulk_run_missing_fq1_exits_cleanly(tmp_path: Path) -> None:
    fq2 = Path("tests/data/bulkdna/L141_CA_R2.fq.gz")
    assert fq2.exists()

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        "L141_CA",
        "--fq1",
        str(tmp_path / "nonexistent_R1.fq.gz"),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--threads",
        "1",
    )
    assert r.returncode == 1, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr
    assert "Forward reads" in r.stderr or "--fq1" in r.stderr


def test_cli_bulk_run_invalid_sample_id_exits_cleanly(tmp_path: Path) -> None:
    fq2 = Path("tests/data/bulkdna/L141_CA_R2.fq.gz")
    assert fq2.exists()

    r = _run(
        "bulk",
        "run",
        "--sample-id",
        "bad/sample",
        "--fq1",
        str(tmp_path / "x.fq.gz"),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--threads",
        "1",
        "--skip-pear",
    )
    assert r.returncode == 1, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr
    assert "sample_id" in r.stderr or "path separators" in r.stderr


def test_cli_bulk_pear_missing_fq2_exits_cleanly(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna/L141_CA_R1.fq.gz")
    assert fq1.exists()

    r = _run(
        "bulk",
        "pear",
        "--sample-id",
        "L141_CA",
        "--fq1",
        str(fq1),
        "--fq2",
        str(tmp_path / "missing_R2.fq.gz"),
        "--output-dir",
        str(tmp_path / "out"),
        "--threads",
        "1",
    )
    assert r.returncode == 1, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr
    assert "Reverse reads" in r.stderr or "--fq2" in r.stderr


def test_bulk_run_skip_pear_accepts_existing_assembled_fastq(tmp_path: Path) -> None:
    fq1 = Path("tests/data/bulkdna/L141_CA_R1.fq.gz")
    fq2 = Path("tests/data/bulkdna/L141_CA_R2.fq.gz")
    assert fq1.exists()
    assert fq2.exists()

    outdir = tmp_path / "out"
    sample_id = "L141_CA"

    pear_run = _run(
        "bulk",
        "pear",
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
    )
    assert pear_run.returncode == 0, f"stdout:\n{pear_run.stdout}\n\nstderr:\n{pear_run.stderr}"

    assembled_fastq = outdir / sample_id / "pear" / "pear.assembled.fastq"
    assert assembled_fastq.exists()

    run = _run(
        "bulk",
        "run",
        "--sample-id",
        sample_id,
        "--output-dir",
        str(outdir),
        "--threads",
        "1",
        "--skip-pear",
        "--assembled-fq",
        str(assembled_fastq),
        "--sample-n",
        "200",
    )
    assert run.returncode == 0, f"stdout:\n{run.stdout}\n\nstderr:\n{run.stderr}"
    assert "Traceback" not in run.stderr
    assert (outdir / sample_id / "reads_1_u_1_l_0.01" / "alleles_by_umis.tsv").exists()


def test_bulk_extract_help_no_longer_mentions_skip_pear() -> None:
    r = _run("bulk", "extract", "--help")
    assert r.returncode == 0
    assert "--skip-pear" not in r.stdout


def test_bulk_finalize_keeps_unannotated_rows_and_logs_counts(tmp_path: Path) -> None:
    combo = ComboPaths(
        dir=tmp_path / "reads_1_u_1_l_0.01",
        denoised_agg_tsv=tmp_path / "reads_1_u_1_l_0.01" / "denoised_agg.tsv",
        denoised_barcodes_tsv=tmp_path / "reads_1_u_1_l_0.01" / "denoised_barcodes.tsv",
        annotated_tsv=tmp_path / "reads_1_u_1_l_0.01" / "annotated.tsv",
        alleles_tsv=tmp_path / "reads_1_u_1_l_0.01" / "alleles_by_umis.tsv",
    )
    combo.ensure_dir()

    denoised_barcodes = pd.DataFrame(
        [
            {"LR": "AAA", "UMIs": 5, "query": "good_query"},
            {"LR": "CCC", "UMIs": 7, "query": "missing_query"},
        ]
    )
    annotated = pd.DataFrame(
        [
            {
                "query": "good_query",
                "query_len": 3,
                "scores": 42.0,
                "mutations": "mut1",
                "aligned_query": "AQ1",
                "aligned_ref": "AR1",
            }
        ]
    )
    denoised_barcodes.to_csv(combo.denoised_barcodes_tsv, sep="\t", index=False)
    annotated.to_csv(combo.annotated_tsv, sep="\t", index=False)

    log_stream = io.StringIO()
    logger = logging.getLogger("test.bulk.finalize")
    logger.handlers = []
    logger.setLevel(logging.INFO)
    logger.propagate = False
    handler = logging.StreamHandler(log_stream)
    logger.addHandler(handler)

    step_finalize(
        denoised_barcodes_tsv=combo.denoised_barcodes_tsv,
        annotated_tsv=combo.annotated_tsv,
        sample_id="L141_CA",
        combo=combo,
        logger=logger,
    )

    final_df = pd.read_csv(combo.alleles_tsv, sep="\t")
    assert len(final_df) == 2
    assert set(final_df["LR"]) == {"AAA", "CCC"}
    aaa = final_df.loc[final_df["LR"] == "AAA"].iloc[0]
    ccc = final_df.loc[final_df["LR"] == "CCC"].iloc[0]
    assert int(aaa["UMIs"]) == 5
    assert aaa["mutation"] == "mut1"
    assert aaa["aligned_query"] == "AQ1"
    assert aaa["aligned_ref"] == "AR1"
    assert int(ccc["UMIs"]) == 7
    assert pd.isna(ccc["mutation"])
    assert pd.isna(ccc["aligned_query"])
    assert pd.isna(ccc["aligned_ref"])
    assert "confidence" not in final_df.columns

    log_text = log_stream.getvalue()
    assert "merged_rows=2" in log_text
    assert "rows_with_aligned_query_ref=1" in log_text
    assert "output_rows=2" in log_text
