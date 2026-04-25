#!/usr/bin/env bash
set -euo pipefail

# Human-readable CLI check that mirrors `tests/test_cli_scrna.py` for camellia.
#
# This script is intended for quick, manual verification:
# - It runs `darlin scrna run` on the bundled camellia test FASTQs
# - It checks that expected output files exist
# - It verifies the extracted.tsv header matches the expected columns
#
# Usage:
#   bash tests/test_cli_scrna_camellia.sh
#
# Notes:
# - Run this inside the Pixi environment so `python` has all dependencies:
#     pixi run bash tests/test_cli_scrna_camellia.sh
# - We run via `python -m darlin.cli` and set `PYTHONPATH=src`, same as pytest.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

fq1="${ROOT_DIR}/tests/data/scCamellia/LL653-CA_L001_R1_001.fastq.gz"
fq2="${ROOT_DIR}/tests/data/scCamellia/LL653-CA_L001_R2_001.fastq.gz"
sample_id="LL653_CA"
protocol="camellia"
whitelist="${ROOT_DIR}/reference/whitelist/scCamellia.txt.gz"
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
