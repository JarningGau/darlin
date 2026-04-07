#!/usr/bin/env bash
set -euo pipefail

# Human-readable CLI check that mirrors `tests/test_cli_bulk_dna.py`.
#
# This script is intended for quick, manual verification:
# - It runs `darlin.cli bulk run` on the bundled test FASTQs
# - It checks that expected output files exist
# - It verifies the alleles TSV header contains required columns
#
# Usage:
#   bash tests/test_cli_bulk_dna.sh
#
# Notes:
# - Run this inside the Pixi environment so `python` has all dependencies:
#     pixi run bash tests/test_cli_bulk_dna.sh
# - We run via `python -m darlin.cli` and set `PYTHONPATH=src`, same as pytest.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

fq1="${ROOT_DIR}/tests/data/bulkdna/L141_CA_R1.fq.gz"
fq2="${ROOT_DIR}/tests/data/bulkdna/L141_CA_R2.fq.gz"
sample_id="L141_CA"

tmp_out="temp"

outdir="${tmp_out}/out"

echo "==> Running bulk DNA pipeline"
echo "    sample_id: ${sample_id}"
echo "    fq1:       ${fq1}"
echo "    fq2:       ${fq2}"
echo "    outdir:    ${outdir}"
echo

python -m darlin.cli bulk run \
  --sample-id "${sample_id}" \
  --fq1 "${fq1}" \
  --fq2 "${fq2}" \
  --output-dir "${outdir}" \
  --threads 1 \
  --keep-pear \
  --sample-n 200

sample_dir="${outdir}/${sample_id}"
combo_dir="${sample_dir}/reads_1_u_1_l_0.01"

expect_file() {
  local p="$1"
  if [[ ! -f "${p}" ]]; then
    echo "ERROR: expected output missing: ${p}" >&2
    exit 1
  fi
}

echo
echo "==> Checking expected outputs exist"
expect_file "${sample_dir}/run.log"
expect_file "${sample_dir}/pear/pear.assembled.fastq"
expect_file "${sample_dir}/extracted.tsv"
expect_file "${sample_dir}/filtered.tsv"
expect_file "${combo_dir}/denoised_barcodes.tsv"
expect_file "${combo_dir}/annotated.tsv"

alleles_tsv="${combo_dir}/alleles_by_umis.tsv"
expect_file "${alleles_tsv}"

echo "==> Checking alleles TSV contains required headers"
ALLELES_TSV="${alleles_tsv}" python - <<'PY'
import csv
import os
import sys

alleles_tsv = os.environ["ALLELES_TSV"]
with open(alleles_tsv, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    if reader.fieldnames is None:
        print("ERROR: TSV has no header row", file=sys.stderr)
        sys.exit(1)
    for name in ("md5", "UMIs"):
        if name not in reader.fieldnames:
            print(f"ERROR: missing required column {name!r} in {reader.fieldnames!r}", file=sys.stderr)
            sys.exit(1)
print("OK: TSV headers include 'md5' and 'UMIs'")
PY

echo
echo "OK: bulk DNA CLI run produced expected outputs"
echo "Output directory was: ${outdir}"
