import os
import subprocess
import sys


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([os.path.abspath("src"), env.get("PYTHONPATH", "")]).strip(os.pathsep)
    return subprocess.run(
        [sys.executable, "-m", "darlin.cli", *args],
        env=env,
        text=True,
        capture_output=True,
    )


def test_top_level_help_mentions_commands() -> None:
    r = _run("--help")
    assert r.returncode == 0
    assert "bulk" in r.stdout
    assert "scrna" in r.stdout


def test_bulk_help_mentions_steps() -> None:
    r = _run("bulk", "--help")
    assert r.returncode == 0
    for step in ["run", "pear", "extract", "filter", "denoise", "annotate"]:
        assert step in r.stdout


def test_bulk_run_help_mentions_no_progress() -> None:
    r = _run("bulk", "run", "--help")
    assert r.returncode == 0
    assert "--no-progress" in r.stdout


def test_scrna_help_mentions_steps() -> None:
    r = _run("scrna", "--help")
    assert r.returncode == 0
    for step in ["run", "extract", "denoise", "qc", "annotate"]:
        assert step in r.stdout


def test_scrna_run_help_mentions_protocol_and_sample_n() -> None:
    r = _run("scrna", "run", "--help")
    assert r.returncode == 0
    assert "--protocol" in r.stdout
    assert "--sample-n" in r.stdout
