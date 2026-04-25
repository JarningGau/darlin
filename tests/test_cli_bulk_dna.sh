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
outdir="${tmp_out}/bulk_dna"

# python -m darlin.cli bulk run \
#   --sample-id "${sample_id}" \
#   --fq1 "${fq1}" \
#   --fq2 "${fq2}" \
#   --output-dir "${outdir}" \
#   --threads 1 \
#   --keep-pear \
#   --locus Col1a1 \
#   --umi-ld 1 2  \
#   --lb-hd-relative 0.01 0.02

python -m darlin.cli bulk run \
  --sample-id "${sample_id}" \
  --fq1 "${fq1}" \
  --fq2 "${fq2}" \
  --output-dir "${outdir}" \
  --threads 1 \
  --keep-pear \
  --locus Col1a1 \
  --umi-ld 1  \
  --lb-hd-relative 0.01 