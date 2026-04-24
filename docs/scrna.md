# scRNA Command Reference

## Overview

The `darlin scrna` command group processes single cell lineage-tracing data (from scrna-seq cDNA amplicon library). The principal entrypoint is `darlin scrna run`, which executes the complete workflow:

`extract -> xx`

## Terms
LB: lineage barcode
CB: cell barcode
UB: UMI
LR: corrected lineage barcode
CR: corrected cell barcode
UR: corrected UMI
HD: hamming distance

## Pipeline
1. extract (LB,CB,UB,reads)
2. correct sequencing error 
  2.1 CB,UB correction (LB,CR,UR,reads)
  2.2 LB correction (LR,CR,UR,reads) (Group by CR,UR,LB_len, HD <= error_rate * LB_len>)
3. 

`darlin scrna run`
