# scRNA Command Reference

## Library Structure
### 10xv3
PE28+350 (PE250)
R1: CB+UB (16+12)
R2: Primer5+DARLIN+Primer3 (350)

### camellia
PE350+28 (PE250)
R1: Primer5+DARLIN+Primer3 (350)
R2: CB+UB (8+8)

Built-in whitelist defaults:

- `10xv3`: `reference/whitelist/10xv3.txt.gz`
- `camellia`: `reference/whitelist/scCamellia.txt.gz`

## Overview

The `darlin scrna` command group processes DARLIN single-cell RNA-seq data. The principal entrypoint is `darlin scrna run`, which executes:

`extract -> denoise -> qc -> annotate`

Most routine analyses should use `run`. Individual subcommands are available for recovery, inspection, and parameter tuning.

## Full-Pipeline Command

```bash
darlin scrna run \
  --sample-id LL837_CA \
  --fq1 tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz \
  --fq2 tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz \
  --protocol 10xv3 \
  --output-dir ./output \
  --sample-n 200
```

## Common Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--sample-id` | required | Sample identifier used for naming the output directory under `--output-dir`. |
| `--output-dir` | `./output` | Base output directory. Results are written under `<output-dir>/<sample-id>/`. |
| `--locus` | `Col1a1` | DARLIN locus name used to load the corresponding `darlin-core` amplicon configuration. |
| `--protocol` | `10xv3` | Single-cell protocol preset. |
| `--whitelist` | protocol default | Override the protocol default whitelist path. |
| `--log-level` | `INFO` | Logging verbosity. Accepted values are `DEBUG`, `INFO`, `WARNING`, and `ERROR`. |
| `--no-progress` | off | Disable tqdm progress bars during extract. |
| `--test` | off | Limit processing to approximately the first 2500 read pairs. |
| `--sample-n` | unset | Limit processing to the first `N` read pairs. |
| `--min-bc-len` | `20` | Minimum lineage-barcode length retained for denoise and annotate. |
| `--fq1` / `--fq2` | required for `run`/`extract` | FASTQs for sequencer R1 and R2. The protocol decides which one is the barcode read. |

## Denoise and QC Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--umi-ld` | `1` | UMI clustering threshold used during per-cell UMI correction. |
| `--lb-error-rate` | `0.01` | Relative lineage-barcode error rate used to derive the Hamming-distance threshold per `(CR, UR, LB_len)` group. |
| `--lb-min-hd` | `1` | Minimum lineage-barcode Hamming-distance threshold. |
| `--major-fraction-threshold-molecule` | `0.8` | Minimum read fraction for the major corrected lineage barcode within a `(CR, UR)` molecule. |
| `--reads-umis-ratio-cutoff` | `1` | Minimum `reads / UMIs` ratio per corrected cell barcode during cell-level QC. |
| `--reads-cutoff` | `1` | Minimum read support per molecule retained after QC. |

## Standard Outputs

The pipeline writes the following standard outputs under `<output-dir>/<sample-id>/`:

- `run.log`
- `extracted.tsv`
- `denoised.tsv`
- `qc.tsv`
- `cell_summary.tsv`
- `annotated.tsv` (columns include `CR`, `LR`, `UR`, `reads`, `mutations`, `aligned_LR`, `aligned_ref`)
- `nUMIs_by_cell_and_lineage.tsv` (columns `n_UMIs`, `CR`, `mutation`, `aligned_LR`, `aligned_ref`)
- `diagnostics/extract_lb_length.png`
- `diagnostics/qc_reads_fraction_hist.png`
- `diagnostics/qc_reads_fraction_scatter.png`
- `diagnostics/qc_reads_vs_umis.png`
- `diagnostics/qc_k_cutoff_curve.png`
- `diagnostics/qc_reads_cutoff_retention.png`
- `diagnostics/qc_n_lr_per_cr_hist.png`

## Step-Wise Execution

Each major stage is also exposed directly:

```bash
darlin scrna extract \
  --sample-id LL837_CA \
  --fq1 tests/data/sc10xv3/LL837-skull-CA_1.fastq.gz \
  --fq2 tests/data/sc10xv3/LL837-skull-CA_2.fastq.gz \
  --output-dir ./output \
  --sample-n 200

darlin scrna denoise \
  --sample-id LL837_CA \
  --output-dir ./output

darlin scrna qc \
  --sample-id LL837_CA \
  --output-dir ./output

darlin scrna annotate \
  --sample-id LL837_CA \
  --output-dir ./output
```

Step commands read the standard upstream output by default. You may override those inputs with:

- `darlin scrna denoise --extracted <path>`
- `darlin scrna qc --denoised <path>`
- `darlin scrna annotate --qc <path>`
