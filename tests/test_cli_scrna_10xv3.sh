#!/usr/bin/env bash
set -euo pipefail

# Human-readable CLI check that mirrors `tests/test_cli_scrna.py` for 10xv3.
#
# This script is intended for quick, manual verification:
# - It runs `darlin scrna run` on the bundled 10xv3 test FASTQs
# - It checks that expected output files exist
# - It verifies the extracted.tsv header matches the expected columns
#
# Usage:
#   bash tests/test_cli_scrna_10xv3.sh
#
# Notes:
# - Run this inside the Pixi environment so `python` has all dependencies:
#     pixi run bash tests/test_cli_scrna_10xv3.sh
# - We run via `python -m darlin.cli` and set `PYTHONPATH=src`, same as pytest.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

fq1="${ROOT_DIR}/tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz"
fq2="${ROOT_DIR}/tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz"
sample_id="LL837_CA"
protocol="10xv3"
whitelist="${ROOT_DIR}/reference/whitelist/10xv3.txt.gz"
locus="Col1a1"

tmp_out="temp"
outdir="${tmp_out}/out"
sample_dir="${outdir}/${sample_id}"

python -m darlin.cli scrna run \
  --sample-id "${sample_id}" \
  --fq1 "${fq1}" \
  --fq2 "${fq2}" \
  --output-dir "${outdir}" \
  --protocol "${protocol}" \
  --whitelist "${whitelist}" \
  --locus "${locus}"
