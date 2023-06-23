# MRaQTL
## Expression regulator activity inference with master regulator (MR) and activity QTL analyses
[![DOI](https://zenodo.org/badge/631020275.svg)](https://zenodo.org/badge/latestdoi/631020275)

This package includes all functions required to apply the master regulator (MR) analysis and activity QTL 
analyses as first demonstrated in [Hoskins et al., 2021](https://doi.org/10.1371/journal.pcbi.1009563) and 
then developed into a user-friendly pipeline for STAR Protocols ([Hoskins et al., 2021](https://doi.org/10.1016/j.xpro.2023.102362)).

The overall pipeline is divided into 4 general steps, 3 of which are in this package.
1. Tissue-specific gene co-expression network inference with ARACNe (not implemented in R).
2. Expression regulators activity inferences and data preparation (implemented in this package).
3. Phenotypic master regulator (MR) inference (implemented in this package).
4. Expression and activity QTL (eQTL and aQTL) analyses (implemented in this package).
    
See the articles linked above and the package's reference manual for details, but here's the TL/DR:
- Identifying genes potentially mediating GWAS signals remains an outstanding challenge
- Master regulators (MRs) integrate genetic/environmental info to establish a cell state
- MR trans-QTL analyses reduce multiple testing burden while enriching for relevant genes
- MRaQTL R package streamlines the approach, empowering post-GWAS hypothesis generation
