# Camellia paired-end fixture

`96c_R1.fq.gz` and `96c_R2.fq.gz` contain the first 10,000 synchronized pairs
from the supplied 89,915-pair sample. Both mates retain their original
300-base sequences and qualities; the `camellia-pe250` preset does not force
250-base trimming.

The full originals are preserved locally in the ignored directory
`jobs/scCamellia-pe250/raw_fastq/`. Their compressed-file SHA-256 checksums are:

- R1: `8d27c30dbcec4fef8f9bc5b3198f222de4a889f0b77a57abfa0f796a16df0752`
- R2: `df2b47051dc8f73b8f7a5e69dfed17780009b14be426c356b2085e2314dc9bbc`

The CLI integration test runs PEAR and all downstream stages on this fixture.
