# Reproduction of: Bulk RNASeq  Retinal Organoids, Scientific Data 2022

This repository contains code used to reproduce figures as presented in: 

> Agarwal, D.; Kuhns, R.; Dimitriou, C.N.; Barlow, E.; Wahlin, K.J.; Enke, R.A. Bulk RNA Sequencing Analysis of Developing Human Induced Pluripotent Cell-Derived Retinal Organoids. *Sci Data *2022*, 9*, 759, https://doi.org/10.1038/s41597-022-01853-x.

----

## Source Data

Raw sorted BAM files were downloaded as provided by the authors in the GitHub repo:
> https://github.com/WahlinLab/Organoid_RNAseq_SciData22/tree/main/BAM_Index

## Software Used

| Software | Version |
| --- | --- | 
| [featureCounts](https://subread.sourceforge.net/) | 2.1.1 | 
| [DESeq2](https://github.com/thelovelab/DESeq2) | 1.52.0 | 
| [ashr](https://github.com/stephens999/ashr) | 2.2-63 | 
| [pheatmap](https://cran.rstudio.com/web/packages/pheatmap/index.html) | 1.0.13 |
| [ggplot2](https://ggplot2.tidyverse.org/) |4.0.3|
| [RColorBrewer](https://cran.rstudio.com/web/packages/RColorBrewer/index.html)|1.1-3|
| [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano) | 1.30.0 |
| [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)|4.20.0|
| [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|3.23.1|
| [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)|1.74.0|
