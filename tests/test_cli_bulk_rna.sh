#!/usr/bin/env bash
set -euo pipefail

# Human-readable CLI check for a bulk RNA sample.
#
# This script is intended for quick, manual verification:
# - It runs `darlin.cli bulk run` on the bundled bulk RNA test FASTQs
#
# Usage:
#   bash tests/test_cli_bulk_rna.sh
#
# Notes:
# - Run this inside the Pixi environment so `python` has all dependencies:
#     pixi run bash tests/test_cli_bulk_rna.sh
# - We run via `python -m darlin.cli` and set `PYTHONPATH=src`, same as pytest.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

fq1="${ROOT_DIR}/tests/data/bulkrna/LL583_RA_1.fastq.gz"
fq2="${ROOT_DIR}/tests/data/bulkrna/LL583_RA_2.fastq.gz"
sample_id="LL583_RA"
locus="Rosa"

tmp_out="temp"
outdir="${tmp_out}/out"

python -m darlin.cli bulk run \
  --sample-id "${sample_id}" \
  --fq1 "${fq1}" \
  --fq2 "${fq2}" \
  --output-dir "${outdir}" \
  --threads 1 \
  --keep-pear \
  --locus "${locus}" \
  --umi-ld 1 \
  --lb-hd-relative 0.01
