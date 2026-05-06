# Bulk Command Reference

## Library structure

### `--protocol pe250` (paired-end overlap assembly)

Typical PE250-style libraries:

- **R1:** UMI — primer3 — DARLIN — primer5  
- **R2:** primer5 — DARLIN — primer3 — UMI  

Reads are assembled with PEAR; extraction uses the assembled FASTQ.

### `--protocol pe85-r350` (short R1 + long R2, no overlap)

PE85+350 libraries do not merge with PEAR. Paired FASTQs are read directly:

- **R1:** UMI at the 5′ end  
- **R2:** primer5 — lineage insert — primer3 (primers matched as forward sequences from the `darlinpy` locus config)

## Overview

The `darlin bulk` command group processes bulk DNA/RNA lineage-tracing data. The principal entrypoint is `darlin bulk run`, which executes:

- **`pe250`:** `pear -> extract -> filter -> denoise -> annotate`
- **`pe85-r350`:** `extract (paired) -> filter -> denoise -> annotate` (PEAR is not run)

Most routine analyses should use `run`. Individual subcommands are available for inspection, parameter tuning, or recovery from an interrupted run.

## Full-Pipeline Command

```bash
darlin bulk run \
  --sample-id <sample> \
  --fq1 <reads_R1.fastq.gz> \
  --fq2 <reads_R2.fastq.gz> \
  [options]
```

### Required arguments

| Argument | Meaning |
|----------|---------|
| `--sample-id` | Sample identifier used for naming the output directory under `--output-dir`. Must be a single path segment (no `/` or `\`). |
| `--fq1` / `--fq2` | Paired FASTQs. Required for `pe250` (unless `--skip-pear`) and required for `pe85-r350`. |

### General input and runtime arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--protocol` | `pe250` | `pe250`: PEAR assembly then extract from assembled reads. `pe85-r350`: paired R1/R2 extraction without PEAR (see Library structure). |
| `--output-dir` | `./output` | Base output directory. Results are written under `<output-dir>/<sample-id>/`. |
| `--locus` | `Col1a1` | DARLIN locus name used to load the corresponding `darlinpy` amplicon configuration. |
| `--pear-path` | `pear` | Path to the PEAR executable (used only when PEAR runs). |
| `--threads` | `8` | Number of threads passed to PEAR when PEAR runs (`pe250` without `--skip-pear`). |
| `--log-level` | `INFO` | Logging verbosity. Accepted values are `DEBUG`, `INFO`, `WARNING`, and `ERROR`. |
| `--no-progress` | off | Disable tqdm progress bars during extract/denoise (cleaner logs in batch or CI). |

### Extraction and Filtering Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--umi-len` | `12` | Number of leading bases interpreted as the UMI during extraction. Reads with `N` in the UMI are excluded. |
| `--min-bc-len` | `20` | Minimum lineage-barcode length retained during filtering and annotation. |

### Denoising Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--reads-cutoff` | `1` | Minimum read support for an aggregated `(lineage barcode, UMI)` pair; enforced when denoising (also used in combo output directory names). |
| `--denoise-iter` | `1` | Number of denoising iterations passed to barcode/UMI correction. |
| `--umi-ld` | `1` | One or more UMI edit-distance thresholds. |
| `--lb-hd-relative` | `0.01` | One or more relative lineage-barcode Hamming-distance thresholds. |

When either `--umi-ld` or `--lb-hd-relative` is provided with multiple values, the pipeline evaluates every parameter combination. Each combination is written to:

`<output-dir>/<sample-id>/reads_<reads-cutoff>_u_<umi-ld>_l_<lb-hd-relative>/`

The `lb-hd-relative` part is formatted with a stable 4-significant-digit general format (`g` conversion), e.g. `0.01` stays `0.01`, `0.0001` becomes `0.0001`.

### PEAR control (`pe250` only)

| Argument | Default | Meaning |
|----------|---------|---------|
| `--skip-pear` | off | Reuse an existing assembled FASTQ at `<output-dir>/<sample-id>/pear/pear.assembled.fastq` instead of rerunning PEAR. Not valid with `--protocol pe85-r350`. |
| `--keep-pear` | off | Preserve the intermediate `pear/` directory after the pipeline completes. By default it is removed. |

### Read subsampling

| Argument | Default | Meaning |
|----------|---------|---------|
| `--test` | off | Limit processing to approximately the first 2500 reads (`pe250`: assembled reads; `pe85-r350`: read pairs scanned). |
| `--sample-n` | unset | Limit processing to the first `N` reads (`pe250`: assembled; `pe85-r350`: paired records). |

## Parameter interactions

- `--sample-n` takes precedence over `--test`. If both are supplied, the explicit `N`-read limit is used.
- `--skip-pear` requires the assembled FASTQ to be present at the expected location unless an individual subcommand accepts an explicit replacement path. It cannot be combined with `--protocol pe85-r350`.
- `--assembled-fq` is only meaningful with `--skip-pear` on `pe250`; it is rejected for `pe85-r350`.
- `--reads-cutoff` is applied in the denoise step (not in `filter`) and appears in combo output-directory names.
- `--umi-ld` and `--lb-hd-relative` define a parameter grid rather than a single joint setting when multiple values are supplied.

## Primary Outputs

| Path | Description |
|------|-------------|
| `<output-dir>/<sample-id>/run.log` | Pipeline log file. New entries are appended (not overwritten) when a subcommand runs again. |
| `<output-dir>/<sample-id>/extracted.tsv` | Extracted lineage barcode and UMI table prior to aggregation. |
| `<output-dir>/<sample-id>/filtered.tsv` | Aggregated `(lineage barcode, UMI)` table after barcode-length filtering; includes all read counts before the denoise-time `--reads-cutoff`. |
| `<output-dir>/<sample-id>/reads_<...>/denoised_agg.tsv` | Denoised barcode/UMI table for one parameter combination. |
| `<output-dir>/<sample-id>/reads_<...>/denoised_barcodes.tsv` | Barcode-level summary with `LR` (lineage barcode), `UMIs`, and `query` (reverse complement) for annotation. |
| `<output-dir>/<sample-id>/reads_<...>/annotated.tsv` | Allele annotation output generated by `darlinpy` (written during `annotate`). |
| `<output-dir>/<sample-id>/reads_<...>/alleles_by_umis.tsv` | Final clone-level output for one parameter combination (written during `annotate`). |
| `<output-dir>/<sample-id>/diagnostics/` | PNG diagnostics from `filter` (pre-denoise curves; uses `--reads-cutoff` for cutoff lines). |
| `<output-dir>/<sample-id>/reads_<...>/` | PNG diagnostics for that parameter combo (pre-denoise, post-denoise, post-annotation), alongside TSV outputs. |

## Step-Wise Execution

Advanced users may execute the workflow one stage at a time.

```bash
darlin bulk pear \
  --sample-id L141_CA \
  --fq1 tests/data/bulkdna/L141_CA_R1.fq.gz \
  --fq2 tests/data/bulkdna/L141_CA_R2.fq.gz \
  --output-dir ./output \
  --threads 1

darlin bulk extract \
  --sample-id L141_CA \
  --output-dir ./output \
  --sample-n 200

darlin bulk filter \
  --sample-id L141_CA \
  --output-dir ./output \
  --min-bc-len 20 \
  --reads-cutoff 1

darlin bulk denoise \
  --sample-id L141_CA \
  --output-dir ./output \
  --reads-cutoff 1 \
  --umi-ld 1 \
  --lb-hd-relative 0.01

darlin bulk annotate \
  --sample-id L141_CA \
  --output-dir ./output \
  --reads-cutoff 1 \
  --umi-ld 1 \
  --lb-hd-relative 0.01 \
  --denoised-barcodes ./output/L141_CA/reads_1_u_1_l_0.01/denoised_barcodes.tsv
```

## Subcommand Notes

- `darlin bulk extract` and `darlin bulk denoise` accept `--no-progress` (same behavior as `bulk run`).
- `darlin bulk pear` assembles paired-end reads and writes `pear/pear.assembled.fastq`.
- `darlin bulk extract` can consume the default assembled FASTQ or an explicit `--assembled-fq` path.
- `darlin bulk filter` can consume the default `extracted.tsv` or an explicit `--extracted` path (length filter and aggregation only). `--reads-cutoff` is not applied at filter time; it is used for diagnostic plot cutoff lines and must match the denoise step you plan to run.
- `darlin bulk denoise` can consume the default `filtered.tsv` or an explicit `--filtered` path.
- `darlin bulk annotate` requires `--denoised-barcodes` and writes both `annotated.tsv` and `alleles_by_umis.tsv` under the combo directory (same naming as `darlin bulk run`). The denoised table should include `query` (as written by `denoise`, or `annotate` can add it from `LR`).
