# Changelog

## [Unreleased]

### Added

### Changed

### Removed

## [0.4.0] - 2026-06-02

**0.4.0** is the repository tree at commit
[`083e4dc`](https://github.com/JarningGau/darlin/commit/083e4dc27b58ab6aeb0daf667cdd6e040286e776)
(`083e4dc27b58ab6aeb0daf667cdd6e040286e776`).

### Added

- Bulk: `darlin bulk extract` accepts `--protocol pe85-r350` with `--fq1`/`--fq2` for paired extraction without PEAR (step-wise workflows aligned with `bulk run`).

### Changed

- Tests: drop redundant scRNA step-by-step and overlapping bulk pe85 integration tests; trim parametrized help/validation cases. Default pytest suite ~50% faster with E2E coverage retained via full `run` tests.
- Bulk: scope `--help` flags per step (`pear`, `extract`, `filter`, `denoise`, `annotate`) so unused common options from `bulk run` are not advertised; mistyped flags are rejected at parse time.
- Bulk/scRNA: tighten `--sample-id` to safe single-segment identifiers (letters, digits, `.`, `_`, `-`; rejects `.`, `..`, and whitespace-only values).
- Bulk: invalid `--locus` and misused `--assembled-fq` (without `--skip-pear` on `pe250`) fail with clean CLI errors instead of tracebacks.
- Bulk/scRNA: reject non-positive numeric CLI thresholds at parse time (`--threads`, `--sample-n`, read cutoffs, UMI length, denoise iteration counts, and related tuning flags).

### Removed

## [0.3.0] - 2026-05-28

**0.3.0** is the repository tree at commit
[`5075168`](https://github.com/JarningGau/darlin/commit/5075168362781cdc2d3c84c33db18fd4548a93e4)
(`5075168362781cdc2d3c84c33db18fd4548a93e4`).

### Added

- scRNA: add denoise-step reads-cutoff diagnostic figure `qc1_reads_cutoff.png`
  (saved under `diagnostic_plots/`, generated before applying the
  `--reads-cutoff-per-molecule` filter).

### Changed

- scRNA: rename QC cutoff flags `--reads-umis-ratio-cutoff` → `--k-cutoff` and
  `--reads-cutoff` → `--reads-cutoff-per-cell`; add `--reads-cutoff-per-molecule`.
- scRNA: shift QC diagnostic plot filename prefixes by one to free the `qc1_`
  slot for the reads-cutoff figure:
  `qc1_pcr_chimera.png` → `qc2_pcr_chimera.png`,
  `qc2_capture_oligo_carryover.png` → `qc3_capture_oligo_carryover.png`,
  `qc2_cells_above_k_cutoff.png` → `qc3_cells_above_k_cutoff.png`,
  `qc3_lineage_barcodes_per_cell.png` → `qc4_lineage_barcodes_per_cell.png`.

### Removed

## [0.2.0] - 2026-05-28

**0.2.0** is the repository tree at commit
[`9d516cd`](https://github.com/JarningGau/darlin/commit/9d516cde3789811d591308867ce6a4ecc6c51340)
(`9d516cde3789811d591308867ce6a4ecc6c51340`).

### Added

- Pixi tasks `test-sc` and `test-bulk` for narrower pytest runs (`pixi.toml`).
- This changelog under `docs/CHANGELOG.md`.

### Changed

- Bump **darlin-core** to **1.1.0**; `reference/allele_bank` materials aligned
  with that release.
- scRNA: QC diagnostic figure filenames and layouts (for example
  `fragment_length_distribution.png`, `qc1_pcr_chimera.png`,
  `qc2_capture_oligo_carryover.png`, `qc2_cells_above_k_cutoff.png`,
  `qc3_lineage_barcodes_per_cell.png`); see `docs/scrna.md`.
- scRNA denoise: molecule counts aggregated after correction using `LR`,
  `CR`, `UR`, and `LB_len` only in `step2_denoised.tsv` (see tests and
  `docs/scrna.md`).
- scRNA annotate output: `nUMIs_by_cell_and_lineage.tsv` renamed to
  `step4_final.tsv` (`paths.final_tsv`).
- scRNA standard outputs: step-prefixed TSV names (`step1_extracted.tsv` through
  `step4_final.tsv`, including `step3_qc_capture_oligo_carryover_data.tsv`) and
  diagnostic figures under `diagnostic_plots/`; see `docs/scrna.md` and
  `ScrnaPaths` in `src/darlin/scrna/paths.py`.
- scRNA pipeline step logs: multi-line summaries with thousands separators,
  inline percentages on extract skip/match lines, and section titles aligned
  with analysis stages (PCR chimera removal, capture oligo carryover removal).

### Removed

- scRNA: standalone `qc_reads_cutoff_retention.png` diagnostic (QC outputs
  consolidated).

## [0.1.0] - 2026-05-13

**0.1.0** is the repository tree at commit
[`971d11c`](https://github.com/JarningGau/darlin/commit/971d11cea8fa45471cdf76d9b2def3fb4cc1176e)
(`971d11cea8fa45471cdf76d9b2def3fb4cc1176e`).

---

When you tag a new version, add a dated `## [x.y.z] - YYYY-MM-DD` section above
**Unreleased**, move finished bullets out of **Unreleased**, and bump
`src/darlin/__init__.py`.
