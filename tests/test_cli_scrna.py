import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


def test_cli_scrna_camellia_pe250_run_produces_outputs(tmp_path: Path) -> None:
    fixture = Path("tests/data/scCamellia-pe250")
    r = _run(
        "scrna", "run", "--sample-id", "96c", "--fq1", str(fixture / "96c_R1.fq.gz"),
        "--fq2", str(fixture / "96c_R2.fq.gz"), "--protocol", "camellia-pe250",
        "--output-dir", str(tmp_path), "--sample-n", "10000", "--threads", "2", "--no-progress",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    sample_dir = tmp_path / "96c"
    extracted = pd.read_csv(sample_dir / "step1_extracted.tsv", sep="\t")
    assert not extracted.empty
    assert extracted.columns.tolist() == ["LB", "CB", "UB", "LB_len"]
    assert extracted["CB"].str.len().eq(8).all()
    assert extracted["UB"].str.len().eq(8).all()
    assert (sample_dir / "pear/pear.assembled.fastq").is_file()
    assert "Input read pairs: 10000" in (sample_dir / "run.log").read_text()
    for filename in ["step2_denoised.tsv", "step3_qc.tsv", "step4_annotated.tsv", "step4_final.tsv"]:
        assert not pd.read_csv(sample_dir / filename, sep="\t").empty


@pytest.mark.parametrize("step", ["run", "extract"])
def test_camellia_pe250_missing_pear_exits_cleanly(tmp_path: Path, step: str) -> None:
    fixture = Path("tests/data/scCamellia-pe250")
    r = _run(
        "scrna", step, "--sample-id", "96c", "--fq1", str(fixture / "96c_R1.fq.gz"),
        "--fq2", str(fixture / "96c_R2.fq.gz"), "--protocol", "camellia-pe250",
        "--pear-path", str(tmp_path / "missing-pear"), "--output-dir", str(tmp_path / "out"),
    )
    assert r.returncode == 1
    assert "PEAR executable not found" in r.stderr
    assert "Traceback" not in r.stderr
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("step", ["run", "extract"])
def test_camellia_pe250_failed_pear_exits_cleanly(tmp_path: Path, step: str) -> None:
    fixture = Path("tests/data/scCamellia-pe250")
    pear = tmp_path / "failed-pear"
    pear.write_text("#!/bin/sh\necho assembly-failed\nexit 3\n")
    pear.chmod(0o755)
    r = _run(
        "scrna", step, "--sample-id", "96c", "--fq1", str(fixture / "96c_R1.fq.gz"),
        "--fq2", str(fixture / "96c_R2.fq.gz"), "--protocol", "camellia-pe250",
        "--pear-path", str(pear), "--output-dir", str(tmp_path / "out"), "--sample-n", "2",
    )
    assert r.returncode == 1
    assert "PEAR assembly failed" in r.stderr
    assert "Traceback" not in r.stderr
    assert "assembly-failed" in (tmp_path / "out/96c/pear/pear.log").read_text()


@pytest.mark.parametrize("sample_args,expected", [([], 2500), (["--sample-n", "7"], 7)])
def test_camellia_pe250_extract_test_limits_input(tmp_path: Path, sample_args, expected) -> None:
    fixture = Path("tests/data/scCamellia-pe250")
    r = _run(
        "scrna", "extract", "--sample-id", "96c", "--fq1", str(fixture / "96c_R1.fq.gz"),
        "--fq2", str(fixture / "96c_R2.fq.gz"), "--protocol", "camellia-pe250",
        "--output-dir", str(tmp_path), "--test", "--threads", "1", "--no-progress", *sample_args,
    )
    assert r.returncode == 0, r.stderr
    assert f"Input read pairs: {expected}" in (tmp_path / "96c/run.log").read_text()


@pytest.mark.parametrize("mismatched", [True, False])
def test_camellia_pe250_invalid_or_unassembled_pairs_exit_cleanly(tmp_path: Path, mismatched: bool) -> None:
    fq1, fq2 = tmp_path / "r1.fq", tmp_path / "r2.fq"
    fq1.write_text(f"@a/1\n{'A' * 150}\n+\n{'I' * 150}\n")
    fq2.write_text(f"@{'b' if mismatched else 'a'}/2\n{'C' * 150}\n+\n{'I' * 150}\n")
    r = _run(
        "scrna", "run", "--sample-id", "empty", "--fq1", str(fq1), "--fq2", str(fq2),
        "--protocol", "camellia-pe250", "--output-dir", str(tmp_path / "out"), "--no-progress",
    )
    assert r.returncode == 1
    assert "Traceback" not in r.stderr
    assert ("ID mismatch" if mismatched else "No lineage barcodes extracted") in r.stderr
    assert not (tmp_path / "out/empty/step4_final.tsv").exists()


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([os.path.abspath("src"), env.get("PYTHONPATH", "")]).strip(os.pathsep)
    return subprocess.run(
        [sys.executable, "-m", "darlin.cli", *args],
        env=env,
        text=True,
        capture_output=True,
    )


def test_cli_scrna_qc_help_mentions_new_cutoffs() -> None:
    r = _run("scrna", "qc", "--help")
    assert r.returncode == 0
    assert "--k-cutoff" in r.stdout
    assert "--reads-cutoff-per-cell" in r.stdout
    assert "--reads-cutoff-per-molecule" not in r.stdout
    assert "--reads-umis-ratio-cutoff" not in r.stdout
    assert "--reads-cutoff " not in r.stdout


def test_scrna_run_missing_fq1_exits_cleanly(tmp_path: Path) -> None:
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    assert fq2.exists()

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        "LL837_CA",
        "--fq1",
        str(tmp_path / "missing_R1.fastq.gz"),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
    )
    assert r.returncode == 1
    assert "Traceback" not in r.stderr
    assert "Forward reads" in r.stderr or "--fq1" in r.stderr


def test_scrna_run_invalid_sample_id_exits_cleanly(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    assert fq1.exists()
    assert fq2.exists()

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        "bad/sample",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--sample-n",
        "10",
    )
    assert r.returncode == 1
    assert "sample_id" in r.stderr


@pytest.mark.parametrize("sample_id", ["."])
def test_scrna_run_rejects_unsafe_sample_id(tmp_path: Path, sample_id: str) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    assert fq1.exists()
    assert fq2.exists()

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        sample_id,
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--sample-n",
        "10",
    )
    assert r.returncode == 1
    assert "Traceback" not in r.stderr
    assert "sample_id" in r.stderr


def test_scrna_run_unknown_protocol_exits_cleanly(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    assert fq1.exists()
    assert fq2.exists()

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        "LL837_CA",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--protocol",
        "unknown_protocol",
        "--sample-n",
        "10",
    )
    assert r.returncode == 1
    assert "unknown protocol" in r.stderr.lower()


def test_scrna_run_missing_whitelist_override_exits_cleanly(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    assert fq1.exists()
    assert fq2.exists()

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        "LL837_CA",
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(tmp_path / "out"),
        "--whitelist",
        str(tmp_path / "missing_whitelist.txt.gz"),
        "--sample-n",
        "10",
    )
    assert r.returncode == 1
    assert "whitelist" in r.stderr.lower()


def test_cli_scrna_denoise_missing_extracted_exits_cleanly(tmp_path: Path) -> None:
    r = _run(
        "scrna",
        "denoise",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
    )
    assert r.returncode == 1
    assert "Run `darlin scrna extract` first" in r.stderr or "--extracted" in r.stderr


def test_cli_scrna_run_produces_outputs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    sample_id = "LL837_CA"
    outdir = tmp_path / "out"

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        sample_id,
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--output-dir",
        str(outdir),
        "--sample-n",
        "200",
        "--no-progress",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    sample_dir = outdir / sample_id
    for path in [
        sample_dir / "run.log",
        sample_dir / "step1_extracted.tsv",
        sample_dir / "step2_denoised.tsv",
        sample_dir / "step3_qc.tsv",
        sample_dir / "step3_qc_capture_oligo_carryover_data.tsv",
        sample_dir / "step4_annotated.tsv",
        sample_dir / "step4_final.tsv",
    ]:
        assert path.exists(), f"Expected output missing: {path}"

    with (sample_dir / "step1_extracted.tsv").open() as f:
        header = f.readline().strip().split("\t")
    assert header == ["LB", "CB", "UB", "LB_len"]

    denoised_df = pd.read_csv(sample_dir / "step2_denoised.tsv", sep="\t")
    for column in ["LR", "CR", "UR", "LB_len", "reads"]:
        assert column in denoised_df.columns

    grouped_df = pd.read_csv(sample_dir / "step4_final.tsv", sep="\t")
    assert list(grouped_df.columns) == [
        "n_UMIs",
        "CR",
        "mutation",
        "aligned_LR",
        "aligned_ref",
    ]


def test_cli_scrna_rejects_non_positive_numeric_args(tmp_path: Path) -> None:
    r = _run(
        "scrna",
        "run",
        "--sample-id",
        "test_sample",
        "--output-dir",
        str(tmp_path / "out"),
        "--fq1",
        "a.fq",
        "--fq2",
        "b.fq",
        "--reads-cutoff-per-cell",
        "0",
    )
    assert r.returncode == 2, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    assert "Traceback" not in r.stderr
    assert "positive" in r.stderr.lower()


def test_cli_scrna_camellia_run_produces_outputs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/scCamellia/LL653-CA_L001_R1_001.fastq.gz")
    fq2 = Path("tests/data/scCamellia/LL653-CA_L001_R2_001.fastq.gz")
    sample_id = "LL653_CA"
    outdir = tmp_path / "out"

    r = _run(
        "scrna",
        "run",
        "--sample-id",
        sample_id,
        "--fq1",
        str(fq1),
        "--fq2",
        str(fq2),
        "--protocol",
        "camellia",
        "--output-dir",
        str(outdir),
        "--sample-n",
        "200",
        "--no-progress",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    sample_dir = outdir / sample_id
    extracted = pd.read_csv(sample_dir / "step1_extracted.tsv", sep="\t")
    assert not extracted.empty
    assert extracted["CB"].str.len().eq(8).all()
    assert extracted["UB"].str.len().eq(8).all()
    assert (sample_dir / "step4_annotated.tsv").exists()
    assert (sample_dir / "step4_final.tsv").exists()
