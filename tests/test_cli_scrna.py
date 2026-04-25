import os
import subprocess
import sys
from pathlib import Path
import pandas as pd


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([os.path.abspath("src"), env.get("PYTHONPATH", "")]).strip(os.pathsep)
    return subprocess.run(
        [sys.executable, "-m", "darlin.cli", *args],
        env=env,
        text=True,
        capture_output=True,
    )


def _run_extract_fixture(tmp_path: Path) -> Path:
    sample_id = "LL837_CA"
    outdir = tmp_path / "out"
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    r = _run(
        "scrna",
        "extract",
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
    return outdir / sample_id


def _run_denoise_fixture(tmp_path: Path) -> Path:
    sample_dir = _run_extract_fixture(tmp_path)
    r = _run(
        "scrna",
        "denoise",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
        "--sample-n",
        "200",
        "--no-progress",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    return sample_dir


def _run_qc_fixture(tmp_path: Path) -> Path:
    sample_dir = _run_denoise_fixture(tmp_path)
    r = _run(
        "scrna",
        "qc",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    return sample_dir


def _run_camellia_extract_fixture(tmp_path: Path) -> Path:
    sample_id = "LL653_CA"
    outdir = tmp_path / "out"
    fq1 = Path("tests/data/scCamellia/LL653-CA_L001_R1_001.fastq.gz")
    fq2 = Path("tests/data/scCamellia/LL653-CA_L001_R2_001.fastq.gz")
    r = _run(
        "scrna",
        "extract",
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
    return outdir / sample_id


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
    assert "sample_id" in r.stderr or "path separators" in r.stderr


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


def test_scrna_run_small_sample_completes_after_validation(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    whitelist = Path("reference/whitelist/10xv3.txt.gz")
    assert fq1.exists()
    assert fq2.exists()
    assert whitelist.exists()

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
        "--sample-n",
        "10",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"
    sample_dir = tmp_path / "out" / "LL837_CA"
    assert (sample_dir / "annotated.tsv").exists()


def test_cli_scrna_extract_produces_outputs(tmp_path: Path) -> None:
    fq1 = Path("tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz")
    fq2 = Path("tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz")
    sample_id = "LL837_CA"
    outdir = tmp_path / "out"

    r = _run(
        "scrna",
        "extract",
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
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    sample_dir = outdir / sample_id
    extracted = sample_dir / "extracted.tsv"
    plot = sample_dir / "diagnostics" / "extract_lb_length.png"
    assert extracted.exists()
    assert plot.exists()

    with extracted.open() as f:
        header = f.readline().strip().split("\t")
    assert header == ["LB", "CB", "UB", "LB_len"]


def test_cli_scrna_camellia_extract_produces_outputs(tmp_path: Path) -> None:
    sample_dir = _run_camellia_extract_fixture(tmp_path)

    extracted = sample_dir / "extracted.tsv"
    plot = sample_dir / "diagnostics" / "extract_lb_length.png"
    assert extracted.exists()
    assert plot.exists()

    df = pd.read_csv(extracted, sep="\t")
    assert list(df.columns) == ["LB", "CB", "UB", "LB_len"]
    assert not df.empty
    assert df["CB"].str.len().eq(8).all()
    assert df["UB"].str.len().eq(8).all()
    assert df["LB_len"].gt(0).any()


def test_cli_scrna_denoise_produces_required_columns(tmp_path: Path) -> None:
    sample_dir = _run_extract_fixture(tmp_path)

    r = _run(
        "scrna",
        "denoise",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
        "--sample-n",
        "200",
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    denoised = sample_dir / "denoised.tsv"
    assert denoised.exists()
    df = pd.read_csv(denoised, sep="\t")
    for column in ["LB", "CB", "UB", "CR", "UR", "LR", "reads", "LB_len"]:
        assert column in df.columns


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


def test_cli_scrna_qc_writes_tables_and_plots(tmp_path: Path) -> None:
    sample_dir = _run_denoise_fixture(tmp_path)

    r = _run(
        "scrna",
        "qc",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    assert (sample_dir / "qc.tsv").exists()
    assert (sample_dir / "cell_summary.tsv").exists()
    for name in [
        "qc_reads_fraction_hist.png",
        "qc_reads_fraction_scatter.png",
        "qc_reads_vs_umis.png",
        "qc_k_cutoff_curve.png",
        "qc_n_lr_per_cr_hist.png",
    ]:
        assert (sample_dir / "diagnostics" / name).exists()


def test_cli_scrna_annotate_produces_required_columns(tmp_path: Path) -> None:
    sample_dir = _run_qc_fixture(tmp_path)

    r = _run(
        "scrna",
        "annotate",
        "--sample-id",
        "LL837_CA",
        "--output-dir",
        str(tmp_path / "out"),
    )
    assert r.returncode == 0, f"stdout:\n{r.stdout}\n\nstderr:\n{r.stderr}"

    annotated = sample_dir / "annotated.tsv"
    assert annotated.exists()
    df = pd.read_csv(annotated, sep="\t")
    for column in ["CR", "LR", "UR", "reads", "mutations", "md5"]:
        assert column in df.columns


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
        sample_dir / "extracted.tsv",
        sample_dir / "denoised.tsv",
        sample_dir / "qc.tsv",
        sample_dir / "cell_summary.tsv",
        sample_dir / "annotated.tsv",
    ]:
        assert path.exists(), f"Expected output missing: {path}"


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
    extracted = pd.read_csv(sample_dir / "extracted.tsv", sep="\t")
    assert not extracted.empty
    assert extracted["CB"].str.len().eq(8).all()
    assert extracted["UB"].str.len().eq(8).all()
    assert (sample_dir / "annotated.tsv").exists()
