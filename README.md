# DARLIN

`DARLIN` is a computational framework for processing lineage-tracing data generated from DARLIN mice, from raw sequencing reads to clone-level inference. The project provides a unified command-line interface for bulk DNA/RNA, single-cell RNA-seq, and single-cell multi-modal assays.

## Installation

```bash
git clone https://github.com/jarninggau/darlin.git
cd darlin
pixi install
pixi run install-cli
```

The workspace is defined in `pixi.toml` and currently targets Python 3.11 on `linux-64`.

## CLI Overview

```bash
darlin --help
```

Primary subcommands:

- `bulk`: recovery of lineage information from bulk DNA/RNA data
- `scrna`: recovery of lineage information from single-cell RNA-seq data

## Quick Start

For bulk datasets, the standard entrypoint is `darlin bulk run`, which executes the complete workflow:

`pear -> extract -> filter -> denoise -> annotate -> finalize`

```bash
darlin bulk run \
  --sample-id L141_CA \
  --fq1 tests/data/bulkdna/L141_CA_R1.fq.gz \
  --fq2 tests/data/bulkdna/L141_CA_R2.fq.gz \
  --output-dir ./output \
  --threads 1

darlin bulk run \
  --sample-id LL583_RA \
  --fq1 tests/data/bulkrna/LL583_RA_1.fastq.gz \
  --fq2 tests/data/bulkrna/LL583_RA_2.fastq.gz \
  --locus Rosa \
  --output-dir ./output \
  --threads 1
```

This command writes sample-specific results under `./output/L141_CA/`.

Detailed bulk documentation, including argument semantics, parameter interactions, output structure, and step-wise execution, is provided in [docs/bulk.md](docs/bulk.md).

## Development Checks

```bash
pixi run smoke
pixi run compile
pixi run test
```

## Status

The project remains under active development. Command-line interfaces and internal APIs should be regarded as provisional.
