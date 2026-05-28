# Changelog

Notable changes to the `darlin` command-line interface, pipeline outputs, and
documentation in this repository.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
The package version in `src/darlin/__init__.py` should match the latest
**released** section below; newer work stays under **Unreleased** until you bump
the version.

**0.1.0** is the repository tree at commit
[`971d11c`](https://github.com/JarningGau/darlin/commit/971d11cea8fa45471cdf76d9b2def3fb4cc1176e)
(`971d11cea8fa45471cdf76d9b2def3fb4cc1176e`).

## [Unreleased]

### Added

### Changed

### Removed

## [0.2.0] - 2026-05-28

**0.2.0** is the repository tree at commit
[`d3aba3a`](https://github.com/JarningGau/darlin/commit/d3aba3a18a80e721f9cd870636a663b88b2b4bfa)
(`d3aba3a18a80e721f9cd870636a663b88b2b4bfa`).

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

First documented release line; points at commit `971d11c` above.

---

When you tag a new version, add a dated `## [x.y.z] - YYYY-MM-DD` section above
**Unreleased**, move finished bullets out of **Unreleased**, and bump
`src/darlin/__init__.py`.
