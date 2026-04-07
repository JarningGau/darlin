# DARLIN

`DARLIN` is a computational framework for processing lineage tracing data from DARLIN mice, from raw sequencing reads to clone-level inference. The package provides a unified CLI interface for bulk DNA/RNA, single-cell (10xv3), and multi-modal (scCamellia) datasets.

## Installation

```bash
git clone https://github.com/jarninggau/darlin.git
cd darlin
pixi install
pixi run install-cli
```

## Quick Start

## CLI Overview

```bash
darlin --help
Commands:
  bulk        Recover lineage information from bulk DNA/RNA data
  scrna       Recover lineage information from single-cell RNA-seq data
  scmulti     Recover lineage information from single-cell multi-modal data
```

### Bulk DNA/RNA

`bulk` is step-oriented: most users run the full pipeline via `run`, while advanced users can execute individual steps.

```bash
# Full pipeline
darlin bulk run --sample-id SAMPLE --fq1 R1.fastq.gz --fq2 R2.fastq.gz --locus Col1a1 --output-dir ./output

# Step-by-step (example)
darlin bulk pear --sample-id SAMPLE --fq1 R1.fastq.gz --fq2 R2.fastq.gz --output-dir ./output
darlin bulk extract --sample-id SAMPLE --output-dir ./output
darlin bulk filter --sample-id SAMPLE --output-dir ./output --reads-cutoff 1 --min-bc-len 20
darlin bulk denoise --sample-id SAMPLE --output-dir ./output --umi-ld 1 --lb-hd-relative 0.01
```

### Development checks

```bash
pixi run smoke
pixi run compile
pixi run test
```

## Pipeline Structure

## Output

## Project Structure

## Status

This project is under active development. Interfaces and APIs may change.

