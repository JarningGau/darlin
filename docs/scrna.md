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
| `--k-cutoff` | `1` | Minimum *k* = `reads / UMIs` ratio per corrected cell barcode during cell-level QC. |
| `--reads-cutoff-per-cell` | `1` | Minimum total reads per corrected cell barcode retained after QC. |
| `--reads-cutoff-per-molecule` | `1` | Minimum reads per molecule retained during denoise. |

## Standard Outputs

The pipeline writes the following standard outputs under `<output-dir>/<sample-id>/`:

- `run.log`
- `step1_extracted.tsv`
- `step2_denoised.tsv`
- `step3_qc.tsv`
- `step3_qc_capture_oligo_carryover_data.tsv`
- `step4_annotated.tsv` (columns include `CR`, `LR`, `UR`, `reads`, `mutations`, `aligned_LR`, `aligned_ref`)
- `step4_final.tsv` (columns `n_UMIs`, `CR`, `mutation`, `aligned_LR`, `aligned_ref`)
- `diagnostic_plots/fragment_length_distribution.png` — histogram of matched read sequence lengths (`LB_len`); vertical dashed line at the unedited lineage-barcode length; empty extract shows "No matched reads".
- `diagnostic_plots/qc1_reads_cutoff.png` — denoise 前分子表：上为 reads cutoff vs 分子数（log x/y），下为 reads cutoff vs 保留 reads 比例；红色虚线为 `--reads-cutoff-per-molecule`。
- `diagnostic_plots/qc2_pcr_chimera.png` — 单一 2×1 图：上面是每个分子 `reads_fraction` 的直方图，红色竖线为 `--major-fraction-threshold-molecule`；下面是 `reads_fraction` vs `reads` 的散点图，y 轴 log，并有相同的阈值线。
- `diagnostic_plots/qc3_capture_oligo_carryover.png` — log–log 散点：每个细胞的 reads vs UMIs（`step3_qc_capture_oligo_carryover_data.tsv`），按 *k* = reads/UMIs 分箱 ≤1、≤5、≤10、>10 着色；虚线红线 slope 1；图例标题为 *k* = Reads/UMIs。
- `diagnostic_plots/qc3_cells_above_k_cutoff.png` — 对 *k* cutoffs 1–19，统计 *k* ≥ cutoff 的 cell 数（y 轴在值较大时使用科学计数法）。
- `diagnostic_plots/qc4_lineage_barcodes_per_cell.png` — 每个 cell (`CR`) 的不同谱系记录数 (`n_LR`) 直方图，整数宽度 bin，y 轴 log。

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
