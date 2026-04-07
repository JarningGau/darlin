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
    assert "scmulti" in r.stdout


def test_bulk_help_mentions_steps() -> None:
    r = _run("bulk", "--help")
    assert r.returncode == 0
    for step in ["run", "pear", "extract", "filter", "denoise", "annotate", "finalize"]:
        assert step in r.stdout

