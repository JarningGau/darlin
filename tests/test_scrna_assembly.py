import logging
import random
from pathlib import Path

import pandas as pd
import pytest

from darlin.scrna.paths import get_scrna_paths
from darlin.scrna.protocols import get_protocol
from darlin.scrna.steps import step_extract


def _fastq(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in records))
    return path


def _extract(tmp_path: Path, fq1: Path, fq2: Path, **kwargs):
    paths = get_scrna_paths(tmp_path / "out", "sample")
    paths.ensure_dirs()
    return step_extract(
        fq1=fq1, fq2=fq2, protocol=get_protocol("camellia-pe250"),
        p5_seq="AGCTACGATCGTACGATGCA", p3_seq="TGCACTAGCTAGCATCGATC",
        unedited_bc_len=360, paths=paths, logger=logging.getLogger("assembly-test"),
        max_reads=None, show_progress=False, **kwargs,
    )


def test_real_pear_recovers_split_primers_and_original_r2_barcode(tmp_path: Path):
    rng = random.Random(20)
    lb = "".join(rng.choices("ACGT", k=360))
    fragment = "AGCTACGATCGTACGATGCA" + lb + "TGCACTAGCTAGCATCGATC"
    reverse = fragment[-234:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    fq1 = _fastq(tmp_path / "r1.fq", [("skip/1 comment", "A" * 250), ("hit/1 comment", fragment[:250])])
    fq2 = _fastq(tmp_path / "r2.fq", [("skip/2 comment", "C" * 250), ("hit/2 comment", "AACCGGTTGGTTAACC" + reverse)])
    result = _extract(tmp_path, fq1, fq2)
    assert result[["LB", "CB", "UB"]].to_dict("records") == [
        {"LB": lb, "CB": "AACCGGTT", "UB": "GGTTAACC"}
    ]
    assert (tmp_path / "out/sample/pear/pear.log").is_file()


def test_assembly_joins_by_id_and_limits_before_pear(tmp_path: Path, monkeypatch, caplog):
    from darlin.scrna import assembly

    fragment = "AGCTACGATCGTACGATGCA" + "ACGT" * 30 + "TGCACTAGCTAGCATCGATC"
    fq1 = _fastq(tmp_path / "r1.fq", [(f"read{i}/1 metadata", fragment[:100]) for i in range(3)])
    fq2 = _fastq(tmp_path / "r2.fq", [
        ("read0/2 metadata", "A" * 16), ("read1/2 metadata", "C" * 8 + "G" * 8),
        ("read2/2 metadata", "T" * 16),
    ])

    def fake_pear(in_fq1, in_fq2, out_prefix, log_file, **kwargs):
        assert Path(in_fq1).read_text().count("\n+") == 2
        assert Path(in_fq2).read_text().count("\n+") == 2
        _fastq(Path(str(out_prefix) + ".assembled.fastq"), [("read1/1 metadata", fragment)])

    monkeypatch.setattr(assembly, "assemble_pe_reads", fake_pear)
    paths = get_scrna_paths(tmp_path / "out", "sample")
    paths.ensure_dirs()
    with caplog.at_level(logging.INFO):
        result = step_extract(
            fq1=fq1, fq2=fq2, protocol=get_protocol("camellia-pe250"),
            p5_seq=fragment[:20], p3_seq=fragment[-20:], unedited_bc_len=120,
            paths=paths, logger=logging.getLogger("assembly-test"),
            max_reads=2, show_progress=False,
        )
    assert result["CB"].tolist() == ["C" * 8]
    assert result["UB"].tolist() == ["G" * 8]
    assert "Input read pairs: 2" in caplog.text
    assert "Assembled read pairs: 1" in caplog.text
    assert "Unassembled read pairs: 1" in caplog.text


@pytest.mark.parametrize("r1_ids,r2_ids,reason", [
    (["a/1"], ["b/2"], "mismatch"),
    (["a/1", "a/1"], ["a/2", "a/2"], "Duplicate"),
    (["a/1"], ["a/2", "b/2"], "different numbers"),
    (["a/1", "b/1"], ["a/2"], "different numbers"),
])
def test_assembly_rejects_invalid_pairing(tmp_path: Path, r1_ids, r2_ids, reason):
    fq1 = _fastq(tmp_path / "r1.fq", [(name, "A" * 150) for name in r1_ids])
    fq2 = _fastq(tmp_path / "r2.fq", [(name, "C" * 150) for name in r2_ids])
    with pytest.raises(ValueError, match=reason):
        _extract(tmp_path, fq1, fq2)


@pytest.mark.parametrize("records,reason", [
    ([("unknown", "ACGT" * 30)], "original pair"),
    ([("a", "ACGT" * 30), ("a", "ACGT" * 30)], "Duplicate assembled"),
])
def test_assembly_rejects_invalid_output_ids(tmp_path: Path, monkeypatch, records, reason):
    from darlin.scrna import assembly

    fq1 = _fastq(tmp_path / "r1.fq", [("a/1", "A" * 150)])
    fq2 = _fastq(tmp_path / "r2.fq", [("a/2", "C" * 150)])

    def fake_pear(in_fq1, in_fq2, out_prefix, log_file, **kwargs):
        _fastq(Path(str(out_prefix) + ".assembled.fastq"), records)

    monkeypatch.setattr(assembly, "assemble_pe_reads", fake_pear)
    with pytest.raises(ValueError, match=reason):
        _extract(tmp_path, fq1, fq2)


def test_zero_assembled_reads_writes_empty_extract(tmp_path: Path):
    fq1 = _fastq(tmp_path / "r1.fq", [("a/1", "A" * 150)])
    fq2 = _fastq(tmp_path / "r2.fq", [("a/2", "C" * 150)])
    result = _extract(tmp_path, fq1, fq2)
    assert result.empty
    extracted = pd.read_csv(tmp_path / "out/sample/step1_extracted.tsv", sep="\t")
    assert extracted.columns.tolist() == ["LB", "CB", "UB", "LB_len"]
