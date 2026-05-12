# Allele Bank

![](figure/allele_bank_experiment.png)

An **allele bank** is a reference collection of observed DARLIN alleles and their associated frequencies. It is used to identify reliable clones by filtering out homoplastic alleles, which can arise independently in multiple cells and therefore confound lineage interpretation.

This allele bank was generated from bulk RNA-seq libraries prepared for each of the three target arrays using 300,000 granulocytes from each of the three DARLIN mouse replicates. Libraries were collected shortly after Dox treatment, at day 3, so allele abundance primarily reflects the frequency at which alleles were generated rather than differential expansion among labeled clones.

## Reproducing the Allele Bank

Preprocess raw FASTQ files into `allele_by_umis_table` outputs:

```bash
./scripts/01.run_DARLIN_bulk_*.py
```

Generate the granulocyte allele bank:

```bash
./scripts/02.generating_allele_bank_Gr.ipynb
```