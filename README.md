# DARLIN

DARLIN is a computational framework for processing lineage-tracing data generated from DARLIN mice, from raw sequencing reads to clone-level inference. The project provides a unified command-line interface for bulk DNA/RNA and single-cell RNA-seq assays.

## Installation

```bash
pixi install
pixi run install-cli
```

The workspace is defined in `pixi.toml` and currently targets Python 3.11 on `linux-64`.

## CLI Overview

```bash
darlin --help
```

Primary subcommands:

| Subcommand | Description | Pipeline |
|------------|-------------|----------|
| `bulk` | Lineage recovery from bulk DNA/RNA data | `pear -> extract -> filter -> denoise -> annotate` |
| `scrna` | Lineage recovery from single-cell RNA-seq data | `extract -> denoise -> qc -> annotate` |

## Quick Start

### Bulk DNA/RNA

```bash
# Bulk DNA (default locus: Col1a1)
darlin bulk run \
  --sample-id L141_CA \
  --fq1 tests/data/bulkdna/L141_CA_R1.fq.gz \
  --fq2 tests/data/bulkdna/L141_CA_R2.fq.gz \
  --output-dir ./output \
  --threads 1

# Bulk RNA (Rosa locus)
darlin bulk run \
  --sample-id LL583_RA \
  --fq1 tests/data/bulkrna/LL583_RA_1.fastq.gz \
  --fq2 tests/data/bulkrna/LL583_RA_2.fastq.gz \
  --locus Rosa \
  --output-dir ./output \
  --threads 1

# Bulk RNA (Tigre locus)
darlin bulk run \
  --sample-id LL583_TA \
  --fq1 tests/data/bulkrna/LL583_TA_1.fastq.gz \
  --fq2 tests/data/bulkrna/LL583_TA_2.fastq.gz \
  --locus Tigre \
  --output-dir ./output \
  --threads 1
```

Results are written under `<output-dir>/<sample-id>/`. See [docs/bulk.md](docs/bulk.md) for the full argument reference, parameter interactions, output layout, and step-wise execution.

### Single-Cell RNA-seq

Supported protocols: `10xv3`, `camellia`.

```bash
# 10x Chromium v3
darlin scrna run \
  --sample-id LL837_CA \
  --fq1 tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz \
  --fq2 tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz \
  --protocol 10xv3 \
  --output-dir ./output \
  --sample-n 200

# Camellia
darlin scrna run \
  --sample-id LL837_CA \
  --fq1 tests/data/sccamellia/LL837-skull-CA_1.fastq.gz \
  --fq2 tests/data/sccamellia/LL837-skull-CA_2.fastq.gz \
  --protocol camellia \
  --output-dir ./output \
  --sample-n 200
```

Results are written under `<output-dir>/<sample-id>/`. See [docs/scrna.md](docs/scrna.md) for the full argument reference, output layout, and step-wise execution.

## Project Layout

```
darlin/
├── src/darlin/
│   ├── cli.py              # Entry point (argparse)
│   ├── commands/            # Subcommand definitions
│   │   ├── bulk.py
│   │   └── scrna.py
│   ├── bulk/                # Bulk DNA/RNA pipeline
│   │   ├── pipeline.py
│   │   ├── steps.py
│   │   ├── denoise.py
│   │   ├── io.py
│   │   ├── matching.py
│   │   ├── pear.py
│   │   ├── plots.py
│   │   └── ...
│   └── scrna/               # scRNA pipeline
│       ├── pipeline.py
│       ├── steps.py
│       ├── protocols.py
│       ├── io.py
│       ├── matching.py
│       ├── plots.py
│       └── ...
├── tests/
├── docs/
├── pixi.toml
└── pyproject.toml
```

## Development

### Pixi Tasks

| Task | Command | Purpose |
|------|---------|---------|
| `smoke` | `pixi run smoke` | CLI help / argparse smoke checks |
| `compile` | `pixi run compile` | `compileall` over `src/` and `tests/` |
| `test` | `pixi run test` | Run pytest |
| `install-cli` | `pixi run install-cli` | Editable pip install (`pip install -e .`) |

### Key Dependencies

| Package | Channel | Purpose |
|---------|---------|---------|
| `pear` | bioconda | Paired-end read assembler (bulk pipeline) |
| `umi_tools` | bioconda | UMI deduplication utilities |
| `darlinpy` | PyPI (git) | DARLIN amplicon annotation engine |
| `biopython` | conda-forge | Sequence I/O |
| `pandas` | conda-forge | Tabular data processing |

All dependencies are managed via Pixi (`pixi.toml` + `pixi.lock`).

## Status

The project remains under active development. Command-line interfaces and internal APIs should be regarded as provisional.
