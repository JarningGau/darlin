# scrna Camellia Protocol Design

Date: 2026-04-25

## Summary

Add `camellia` as a built-in `darlin scrna` protocol. Camellia uses the same downstream `denoise -> qc -> annotate` behavior as `10xv3`; only extraction defaults and read layout differ.

## Protocol Contract

- Protocol name: `camellia`
- Cell barcode length: `8`
- UMI length: `8`
- Default whitelist: `reference/whitelist/scCamellia.txt.gz`
- Barcode read: `fq2`
- DARLIN amplicon read: `fq1`

## Implementation Notes

- Extend `ScrnaProtocol` so read roles are declarative protocol metadata rather than special-cased branches.
- Update extraction to read `CB` and `UB` from the configured barcode read and perform primer matching on the configured DARLIN read.
- Keep `denoise`, `qc`, and `annotate` protocol-agnostic and unchanged.

## Verification

- Add CLI tests for `scrna extract --protocol camellia` using `tests/data/scCamellia`.
- Add a small end-to-end `scrna run --protocol camellia` test that verifies extraction lengths and final outputs.
