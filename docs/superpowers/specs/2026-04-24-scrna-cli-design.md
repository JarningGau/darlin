# scrna CLI Design

Date: 2026-04-24

## Summary

Implement a production `darlin scrna` command group aligned with the existing `darlin bulk` CLI style. The command group will expose the major single-cell RNA-seq processing stages as subcommands:

- `run`
- `extract`
- `denoise`
- `qc`
- `annotate`

The implementation will productize the workflow already proven in `notebook/darlin-scrna.py` and preserve its current behavior for the `10xv3` protocol. The design must also leave a clear extension point for future protocols without rewriting the main pipeline.

## Goals

- Match the `bulk` command style and recovery model.
- Keep the command layer thin and move processing logic into `src/darlin/scrna/`.
- Support step-wise execution and full-pipeline execution.
- Preserve the notebook workflow:
  `extract -> denoise -> qc -> annotate`
- Emit the QC plots currently produced in the notebook.
- Keep the first implementation easy to test with existing bundled `sc10xv3` data and low-friction assertions.
- Add protocol-aware configuration so `10xv3` is the default, not a hard-coded special case spread across the pipeline.

## Non-Goals

- Supporting multiple protocols in this change beyond a clean `10xv3`-first extension point.
- Adding new biological filtering rules beyond what the notebook already implements.
- Redesigning `bulk` internals or forcing code sharing that would distort either workflow.

## Command Surface

`darlin scrna` will be a command group parallel to `darlin bulk`.

### Subcommands

- `darlin scrna run`
  Runs the full pipeline in order.
- `darlin scrna extract`
  Extracts lineage barcode, cell barcode, and UMI information from paired FASTQs.
- `darlin scrna denoise`
  Collapses reads and performs CB, UMI, and lineage-barcode correction.
- `darlin scrna qc`
  Applies molecular-level and cell-level QC filters and writes QC diagnostics.
- `darlin scrna annotate`
  Annotates corrected lineage barcodes with `darlinpy`.

### General CLI Principles

- Required and optional arguments should read like `bulk`.
- Every step reads a standard upstream output by default and may accept an explicit replacement input path.
- User-facing validation errors should exit cleanly without stack traces.
- `run` should be the routine entrypoint; step commands are for recovery, inspection, and parameter tuning.

## Module Layout

Add a new package under `src/darlin/scrna/` with focused modules:

- `protocols.py`
  Protocol definitions and lookup.
- `paths.py`
  Sample output path resolution.
- `validate.py`
  Input validation helpers.
- `logging.py`
  Log setup aligned with `bulk`.
- `io.py`
  FASTQ and whitelist helpers.
- `matching.py`
  Primer matching and small sequence utilities.
- `steps.py`
  Step implementations for `extract`, `denoise`, `qc`, and `annotate`.
- `plots.py`
  QC figure generation and output.
- `pipeline.py`
  Full-pipeline orchestration for `run`.

`src/darlin/commands/scrna.py` should remain thin and limited to parser wiring, argument normalization, path resolution, input validation, and dispatch.

## Protocol Design

Protocol must be a first-class input instead of an implementation detail.

### Initial CLI Contract

- `--protocol`
  Default: `10xv3`
- `--whitelist`
  Optional override for the protocol default whitelist path

### Initial Protocol Model

The first version only needs enough structure to avoid hard-coding `10xv3` logic in each step. A protocol definition should include at least:

- `name`
- `cb_len`
- `umi_len`
- `default_whitelist_path`

If later protocols need different R1 parsing behavior, the model can grow to include a protocol-specific parser function or parser metadata. The main pipeline should consume the protocol definition rather than branch on protocol names inline.

### Initial Built-In Protocol

`10xv3`:

- `cb_len = 16`
- `umi_len = 12`
- default whitelist path: `reference/whitelist/10xv3.txt.gz`

## Data Flow

The pipeline follows the notebook implementation:

`extract -> denoise -> qc -> annotate`

### 1. Extract

Inputs:

- `fq1`
- `fq2`
- `locus`
- `protocol`

Behavior:

- Load primer definitions from `darlinpy` using `locus`.
- Parse R1 using the selected protocol to split `CB` and `UB`.
- Search R2 for one `p5` and one `p3` match using the notebook logic.
- Extract `LB` from the sequence between the primer matches.
- Skip records with ambiguous bases in `CB` or `UB`.
- Record extraction summary statistics in the log.

Primary output:

- `<output-dir>/<sample-id>/extracted.tsv`

Columns:

- `LB`
- `CB`
- `UB`
- `LB_len`

Diagnostics:

- lineage-barcode length histogram PNG in `diagnostics/`

### 2. Denoise

Input:

- default: `<sample>/extracted.tsv`
- override: `--extracted`

Behavior:

- Collapse duplicate `(LB, CB, UB)` observations to a `reads` count.
- Apply minimum lineage-barcode length filtering.
- Correct `CB -> CR` against the whitelist using the notebook’s HD1-neighbor semantics.
- Correct `UB -> UR` per corrected cell barcode using `umi_tools.UMIClusterer`.
- Correct `LB -> LR` within `(CR, UR, LB_len)` groups using the notebook error-rate rule:
  `max(round(error_rate * LB_len), min_hd)`
- Re-aggregate after each correction stage as needed.
- Log summary counts after major denoising phases.

Primary output:

- `<output-dir>/<sample-id>/denoised.tsv`

Columns:

- `LB`
- `CB`
- `UB`
- `CR`
- `UR`
- `LR`
- `reads`
- `LB_len`

The denoised table intentionally retains pre-correction and post-correction identifiers to make downstream debugging easier.

### 3. QC

Input:

- default: `<sample>/denoised.tsv`
- override: `--denoised`

Behavior:

- Molecular QC:
  compute total reads per `(CR, UR)`, derive `reads_fraction`, and keep only major `LR` assignments meeting `major_fraction_threshold_molecule`.
- Cellular QC:
  summarize each `CR` with `n_reads`, `n_UR`, and `k = n_reads / n_UR`.
- Filter cells using the notebook thresholds:
  `k >= reads_umis_ratio_cutoff`
  and molecule rows with `reads >= reads_cutoff`
- Compute `n_LR` per final `CR` for plotting.
- Log the number of reads removed by molecular and cellular QC.

Primary outputs:

- `<output-dir>/<sample-id>/qc.tsv`
- `<output-dir>/<sample-id>/cell_summary.tsv`

`qc.tsv` columns:

- `CR`
- `UR`
- `LR`
- `reads`
- `k`
- `n_LR`

Additional columns may be retained if convenient for debugging, but these fields must exist.

Diagnostics:

- histogram of `reads_fraction`
- scatter of `reads_fraction` vs `reads`
- scatter of `n_reads` vs `n_UR`
- line plot of number of `CR`s above each `k` cutoff
- histogram of `n_LR` per `CR`

All QC PNGs are written under `<sample>/diagnostics/`.

### 4. Annotate

Input:

- default: `<sample>/qc.tsv`
- override: `--qc`

Behavior:

- Deduplicate `LR` values from QC output.
- Call `darlinpy` annotation on unique corrected lineage barcodes.
- Merge annotation results back onto molecule rows.
- Compute `md5` from `aligned_query + aligned_ref` using the notebook rule.

Primary output:

- `<output-dir>/<sample-id>/annotated.tsv`

Columns:

- `CR`
- `LR`
- `UR`
- `reads`
- `mutations`
- `md5`

Additional annotation columns may be preserved if helpful, but the final table must contain the columns above.

## Output Layout

The directory layout should follow the existing `bulk` conventions:

- `<output-dir>/<sample-id>/run.log`
- `<output-dir>/<sample-id>/extracted.tsv`
- `<output-dir>/<sample-id>/denoised.tsv`
- `<output-dir>/<sample-id>/qc.tsv`
- `<output-dir>/<sample-id>/cell_summary.tsv`
- `<output-dir>/<sample-id>/annotated.tsv`
- `<output-dir>/<sample-id>/diagnostics/*.png`

Each step should create required directories before writing outputs.

## Parameters

Expose only parameters already justified by the notebook and current workflow.

### Common Parameters

- `--sample-id`
- `--output-dir` default `./output`
- `--locus` default `Col1a1`
- `--log-level` with the same choices as `bulk`
- `--no-progress`

### Input Parameters

- `run` and `extract`:
  `--fq1`, `--fq2`
- `denoise`:
  `--extracted`
- `qc`:
  `--denoised`
- `annotate`:
  `--qc`

### Protocol Parameters

- `--protocol` default `10xv3`
- `--whitelist` optional override

### Denoise and QC Parameters

- `--min-bc-len`
- `--umi-ld`
- `--lb-error-rate`
- `--lb-min-hd`
- `--major-fraction-threshold-molecule`
- `--reads-umis-ratio-cutoff`
- `--reads-cutoff`

### Low-Friction Development Parameters

- `--test`
  Limit processing to approximately the first 2500 read pairs.
- `--sample-n`
  Limit processing to the first `N` read pairs.

As in `bulk`, `--sample-n` takes precedence over `--test`.

## Logging and Errors

- Reuse the `bulk` logging model: write to `<sample>/run.log` and append across reruns.
- Each step should log a concise summary with counts relevant to that stage.
- Validation errors should be converted to single-line or short multi-line user-facing messages on `stderr`.
- Missing upstream outputs should guide the user to either run the previous step or pass the explicit override flag.

## Testing Strategy

The implementation should optimize for low-friction tests rather than fragile exact-value assertions.

### 1. Smoke Tests

Extend CLI smoke coverage to verify:

- `darlin --help` still lists `scrna`
- `darlin scrna --help` lists all major steps
- step help includes the expected key options such as `--protocol` and `--sample-n`

### 2. Step-Wise CLI Tests

Use bundled `tests/data/sc10xv3/*` inputs with small `--sample-n` values.

Test cases should cover:

- `extract` writes `extracted.tsv` and at least one PNG
- `denoise` writes `denoised.tsv` with required columns
- `qc` writes `qc.tsv`, `cell_summary.tsv`, and all expected QC plots
- `annotate` writes `annotated.tsv` with `mutations` and `md5`

Assertions should focus on:

- process exits cleanly
- files exist
- required columns exist
- logs contain expected stage summaries

### 3. Minimal End-to-End Run Test

Add one full `darlin scrna run` test using bundled `sc10xv3` data and small `--sample-n`.

Assertions should verify:

- full pipeline exits with return code `0`
- all primary outputs exist
- diagnostic PNGs exist
- `run.log` contains entries for extract, denoise, qc, annotate, and final timing summary

### 4. Validation Tests

Add a few clean failure-path tests similar to `bulk`:

- missing FASTQ
- invalid `sample_id`
- missing upstream intermediate file for a step command
- unknown protocol
- missing whitelist override path

## Implementation Notes

- Follow the notebook behavior first; avoid opportunistic biological-method changes in this pass.
- Prefer small helper functions with explicit inputs over hidden module state.
- Keep plotting isolated from filtering logic so tests can verify the data outputs separately from image generation.
- Avoid over-generalizing protocol support in the first pass; the extension point matters more than speculative abstraction.

## Acceptance Criteria

- `darlin scrna` is a working command group with `run`, `extract`, `denoise`, `qc`, and `annotate`.
- The default workflow on `10xv3` reproduces the notebook’s major processing stages.
- Each step can be run independently using standard default inputs or explicit override paths.
- QC plots are written as PNG files.
- The design leaves a clean path to add more protocols by extending protocol definitions rather than rewriting step logic.
- Tests cover help text, step-wise execution, one full pipeline run, and common validation failures.
