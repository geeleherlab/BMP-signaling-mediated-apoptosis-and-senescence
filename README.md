# Introduction

Here are the R scripts used for high-throughput ChIP-seq and RNA-seq data analysis in the research study titled 
'Bone Morphogenetic Protein (BMP) Signaling in Determining Neuroblastoma Cell Fate and Sensitivity to Retinoic Acid'.


# Required packages and setup environment 

```sh
R version 4.1.2
system x86_64

library required:
"ggplot2"
"tidyverse"
"DESeq2"
"clusterProfiler"
"enrichplot"
"ggvenn"
"cowplot"
```

# Scripts for ChIP-seq analysis  

`ChIP.genome.anno.r` the script use ChIP-seq peaks, the narrowpeak files that called from mac2, as input files, and provide annotations of distribution of transciption factor binding loci relative to transcription start site.  

`ChIP.GSEApathway.venn.r` The script utilizes the combined peaks from various transcription factors, such as RARA and SMAD4, to calculate the overlapping number of peaks among different transcription factors. Additionally, the script analyzes the differential ChIP signals for each transcription factor under specific conditions. It also performs GSEA annotation based on peak intensity for each ChIP-seq samples.  

`ChIP.overlapTF.GSEA.r` The script conducts GSEA and hypergeometric-based pathway enrichment analysis on overlapping peaks among various transcription factors in a specific cell type. An example includes analyzing the binding peaks between transcription factors RARA and SMAD9 in the CHP-134 cell line.

`ChIP.peak.correlation.r` The script calculates the correlation matrix for the binding intensity of peaks between each pair of samples, considering various treatment conditions in a specific cell type.

# Scripts for RNA-seq analysis

# Scripts for integrative analysis between ChIP-seq and RNA-seq data

# Scripts for single cell RNA-seq data comparison
