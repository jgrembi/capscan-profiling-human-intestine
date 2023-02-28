# capscan-profiling-human-intestine
Profiling of the human intestinal microbiome and bile acids with CapScan sampling device
[Zenodo DOI](10.5281/zenodo.7683655)

## Description
This repository includes data and replication files to reproduce the analyses in the manuscript:

Shalon, D*, Culver, RN*, Grembi, JA*, Folz, J*, Treit, PV, Shi, H, Rosenberger, FA, Dethlefsen, L, Meng, X, Yaffe, E, Aranda-Diaz, A, Geyer, PE, Mueller-Reif, JB, Spencer, S, Patterson, AD, Triadafilopoulos, G, Holmes, SP, Mann, M, Fiehn, O, Relman, DA, Huang, KC. Profiling the human intestinal environment under physiologic conditions.  _Nature_. (2023)


The configuration file `0-config.R` lists all R packages that need to be installed prior to running the code. 

The `data` directory includes all datasets needed to reproduce the analysis plus a data dictionary (`0-data-dictionary-samples.xlsx`), which decribes the metadata variables present in the files `sample_data.RDS`, `phyloseq_alphaDiv.rds`, and `phyloseq_bilesalt.rds`. The latter two are phyloseq objects with data in each of the object slots (sample data, ASV table, taxonomy table, and phylogenetic tree). 

Raw data is available from the following public repositories:
 - Sequencing data (both 16S and shotgun metagenomic): NCBI Sequencing Read Archive under [BioProject PRJNA822660](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA822660)
 - Metabolomics: Metabolomics Workbench under project numbers [ST002073](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001315) and [ST002075](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001315)
 - Proteomics: ProteomeXchange Consortium via the PRIDE partner repository with dataset identifier [PXD038906](https://www.ebi.ac.uk/pride/archive/projects/PXD038906)


Questions can be directed to Jess Grembi (jgrembi@stanford.edu) or KC Huang (kchuang@stanford.edu).
