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
```

# Scripts for ChIP-seq analysis
`ChIP.genome.anno.r` the script use ChIP-seq peaks, the narrowpeak files that called from mac2, as input files, and provide annotations of distribution of transciption factor binding loci relative to transcription start site.


# Scripts for RNA-seq analysis

# Scripts for integrative analysis between ChIP-seq and RNA-seq data

# Scripts for single cell RNA-seq data comparison
