# scRNA Command Reference

## Library Structure
### 10xv3
PE28+350
R1: CB+UB (16+12)
R2: Primer5+DARLIN+Primer3 (350)

### camellia
PE350+28
R1: Primer5+DARLIN+Primer3 (350)
R2: CB+UB (8+8)

### camellia-pe250

Overlapping paired-end Camellia reads (PE250 or longer, including PE300).
R1 begins with the DARLIN amplicon; R2 begins with CB+UB (8+8) followed by
the opposite end of the amplicon. PEAR assembles the pair before extracting
the lineage barcode. CB and UB are always taken from the original R2,
associated with the assembled read by ID. Unassembled and discarded pairs
do not enter downstream analysis. No fixed-length trimming is applied.

Built-in whitelist defaults:

- `10xv3`: `reference/whitelist/10xv3.txt.gz`
- `camellia`: `reference/whitelist/scCamellia.txt.gz`
- `camellia-pe250`: `reference/whitelist/scCamellia.txt.gz`

## Overview

The `darlin scrna` command group processes DARLIN single-cell RNA-seq data. The principal entrypoint is `darlin scrna run`, which executes:

`extract -> denoise -> qc -> annotate`

For `camellia-pe250`, both `run` and `extract` execute PEAR first:
`pear -> extract -> denoise -> qc -> annotate`.

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

For PEAR-first Camellia processing:

```bash
darlin scrna run \
  --sample-id 96c \
  --fq1 tests/data/scCamellia-pe250/96c_R1.fq.gz \
  --fq2 tests/data/scCamellia-pe250/96c_R2.fq.gz \
  --protocol camellia-pe250 \
  --locus Col1a1 \
  --output-dir ./output \
  --sample-n 10000 \
  --threads 8
```

## Common Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `--sample-id` | required | Sample identifier used for naming the output directory under `--output-dir`. |
| `--output-dir` | `./output` | Base output directory. Results are written under `<output-dir>/<sample-id>/`. |
| `--locus` | `Col1a1` | DARLIN locus name used to load the corresponding `darlin-core` amplicon configuration. |
| `--protocol` | `10xv3` | Single-cell protocol preset: `10xv3`, `camellia`, or `camellia-pe250`. |
| `--whitelist` | protocol default | Override the protocol default whitelist path. |
| `--log-level` | `INFO` | Logging verbosity. Accepted values are `DEBUG`, `INFO`, `WARNING`, and `ERROR`. |
| `--no-progress` | off | Disable tqdm progress bars during extract. |
| `--test` | off | Limit processing to approximately the first 2500 read pairs. |
| `--sample-n` | unset | Limit processing to the first `N` read pairs. |
| `--min-bc-len` | `20` | Minimum lineage-barcode length retained for denoise and annotate. |
| `--fq1` / `--fq2` | required for `run`/`extract` | FASTQs for sequencer R1 and R2. The protocol decides which one is the barcode read. |
| `--pear-path` | `pear` | PEAR executable for `camellia-pe250`; available on `run` and `extract`. |
| `--threads` | `8` | Positive number of PEAR assembly threads; available on `run` and `extract`. |

For `camellia-pe250`, `--test` and `--sample-n` limit input pairs **before**
assembly. `--sample-n` overrides `--test`. PEAR uses its default assembly
parameters; `--pear-path` and `--threads` are unused for other presets.
Read IDs must match across mates and be unique within the selected input.
The first whitespace-delimited ID token is used, without a terminal `/1` or
`/2` mate suffix. Invalid pairing and PEAR failures produce a nonzero exit
status with an explanatory error. An empty assembly writes a header-only
extract table; `run` stops with an error if extraction yields no lineage
barcodes.

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
- `diagnostic_plots/qc1_reads_cutoff.png` — denoise-stage molecule table before filtering: top panel shows reads cutoff vs number of molecules (log–log axes); bottom panel shows reads cutoff vs fraction of reads retained; the red dashed line marks `--reads-cutoff-per-molecule`.
- `diagnostic_plots/qc2_pcr_chimera.png` — single 2×1 figure: the top panel is a histogram of per-molecule `reads_fraction` with a red vertical line at `--major-fraction-threshold-molecule`; the bottom panel is a scatter plot of `reads_fraction` vs `reads` with a log-scaled y axis and the same threshold line.
- `diagnostic_plots/qc3_capture_oligo_carryover.png` — log–log scatter plot of per-cell reads vs UMIs (from `step3_qc_capture_oligo_carryover_data.tsv`), with points colored by *k* = reads/UMIs bins (≤1, ≤5, ≤10, >10); a dashed red line with slope 1; legend title *k* = Reads/UMIs.
- `diagnostic_plots/qc3_cells_above_k_cutoff.png` — for *k* cutoffs 1–19, counts the number of cells with *k* ≥ cutoff (y axis uses scientific notation for large values).
- `diagnostic_plots/qc4_lineage_barcodes_per_cell.png` — histogram of the number of distinct lineage barcodes (`n_LR`) per cell (`CR`), with integer-width bins and a log-scaled y axis.

For `camellia-pe250`, PEAR products are also retained in `pear/`:
`pear.log`, `pear.assembled.fastq`, `pear.unassembled.forward.fastq`,
`pear.unassembled.reverse.fastq`, and `pear.discarded.fastq`.
Temporary normalized input copies are removed after assembly. `run.log`
reports input and assembled pair counts, and the number of pairs not
assembled (including discarded pairs). Extraction counts refer to assembled
reads. Repeating `run` or `extract` runs PEAR again rather than reusing a prior
assembly. Existing TSV filenames and columns remain the same.

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
