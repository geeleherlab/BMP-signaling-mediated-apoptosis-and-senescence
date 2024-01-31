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
# input dataset
The input files required runing scripts below can be asscess here https://osf.io/tx7je/files/osfstorage.

# Scripts

`find.ref.peak.pl` The script processes narrowPeak files produced by MACS2 for each sample and creates a pipeline using BEDTools to calculate the consistent peak coordinates across different samples.  
 
`ChIP.genome.anno.r` the script use ChIP-seq peaks, the narrowpeak files that called from mac2, as input files, and provide annotations of distribution of transciption factor binding loci relative to transcription start site.  

`ChIP.GSEApathway.venn.r` The script utilizes the combined peaks from various transcription factors, such as RARA and SMAD4, to calculate the overlapping number of peaks among different transcription factors. Additionally, the script analyzes the differential ChIP signals for each transcription factor under specific conditions. It also performs GSEA annotation based on peak intensity for each ChIP-seq samples.  

`ChIP.overlapTF.GSEA.r` The script conducts GSEA and hypergeometric-based pathway enrichment analysis on overlapping peaks among various transcription factors in a specific cell type. An example includes analyzing the binding peaks between transcription factors RARA and SMAD9 in the CHP-134 cell line.

`ChIP.peak.correlation.r` The script calculates the correlation matrix for the binding intensity of peaks between each pair of samples, considering various treatment conditions in a specific cell type.

`RNA.pathwayanno.cellline.r` The script executes GSEA and hypergeometric-based pathway enrichment analysis on differentially expressed genes between two samples under comparison, utilizing comprehensive annotation databases such as KEGG, GO, and Bioplanet. Additionally, the script performs PCA analysis on various samples based on the gene expression read counts of all genes.  

`RNA.tumor.cell.compare.r` The script reads gene expression read counts for both cell line and tumor samples. It then employs DEseq2 to calculate fold change, p-value, and differentially expressed genes between these samples. Additionally, the script compares the fold changes observed in cell lines with those between tumor and cell line samples. It also performs GSEA pathway annotation based on the fold change data obtained from comparing cell line and tumor samples.

`Integrative.RNA.ChIP.r` The script perform integrative analysis of RNA-seq and ChIP-seq data for the CHP-134 cell line, it calculates the proportion of differentially expressed genes between treatment and control groups that are also marked by relevant ChIP binding peaks. This analysis includes creating a portion plot and conducting statistical tests to determine the significance of the overlap. Additionally, it calculates the number of differentially expressed genes within each functional group of interest.

`Singlecell.heatmap.r` The script reads the expression levels of target genes from neuroblastoma programs (single cell clusters) identified in single cell RNA-seq data and generates a heatmap to display the expression levels of selected genes across five different single cell RNA-seq datasets.

`Singlecell.boxplot.r` The script analyzes the expression levels of targeted genes within neuroblastoma programs across five distinct single cell RNA-seq datasets, and generates a boxplot to compare the expression levels of these selected genes across the various datasets.


