# capscan-profiling-human-intestine
Profiling of the human intestinal microbiome and bile acids with CapScan sampling device

## Description
This repository includes data and replication files to reproduce the analyses in the manuscript:

Shalon, D*, Culver, RN*, Grembi, JA*, Folz, J*, Treit, PV, Shi, H, Rosenberger, FA, Dethlefsen, L, Meng, X, Yaffe, E, Aranda-Diaz, A, Geyer, PE, Mueller-Reif, JB, Spencer, S, Patterson, AD, Triadafilopoulos, G, Holmes, SP, Mann, M, Fiehn, O, Relman, DA, Huang, KC. _Profiling the human intestinal environment under physiologic conditions._  Nature. (2023)


The `data` directory includes all datasets needed to reproduce the analysis plus a data dictionary (`0-data-dictionary-samples.xlsx`), which decribes the metadata variables present in the files `sample_data.RDS`, `phyloseq_alphaDiv.rds`, and `phyloseq_bilesalt.rds`. The latter two are phyloseq objects with data in each of the object slots (sample data, ASV table, taxonomy table, and phylogenetic tree). The configuration file `0-config.R` lists all R packages that need to be installed prior to running the code. 

Raw data can be found at the following public repositiries:
 - Sequencing data (both 16S and shotgun metagenomic): NCBI Sequencing Read Archive under BioProject PRJNA822660
 - Metabolomics: Metabolomics Workbench under project numbers ST002073 and ST002075
 - Proteomics: ProteomeXchange Consortium via the PRIDE partner repository with dataset identifier PXD038906


Questions can be directed to Jess Grembi (jgrembi@stanford.edu) or KC Huang (kchuang@stanford.edu).
