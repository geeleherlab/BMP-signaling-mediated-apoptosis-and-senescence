setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/1_ChIP/")
#setwd( "/Volumes/groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/chip_consensus_ref/" )
library("ggplot2")
library(tidyverse)
library("clusterProfiler")
library("ChIPseeker")
library("ChIPpeakAnno")
library(cowplot)
#################################################################################CHP-134 SMAD9
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

files=Sys.glob("./data/1_CHP134/peaks/SMAD9/*.narrowPeak")
files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
peaks <- lapply(files2,readPeakFile)
names(peaks)=c("0_1","0_2","1D_2","1_1","3_1","3_2","6_1","6_2")
# peak
###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############b
p1=plotAnnoBar(peakAnnoList)
#p6=upsetplot(peakAnno)

##########################distribution of tf binding loci relative to TSS
p2=plotDistToTSS(peakAnnoList,
                 title="Distribution of transcription factor-binding loci\nrelative to TSS")


pdf("./CHP134.SMAD9.distribution.pdf", width=12, height=5)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()


###############################################################################CHP-134 RARA

files=Sys.glob("./data/1_CHP134/peaks/RARA/Min*.macs2.filter.narrowPeak")
files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}

peaks <- lapply(files2,readPeakFile)
names(peaks)=c("DMSO_1","DMSO_2","RA_1","RA_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)
#p6=upsetplot(peakAnno)

##########################distribution of tf binding loci relative to TSS
p2=plotDistToTSS(peakAnnoList,
                  title="Distribution of transcription factor-binding loci\nrelative to TSS")


pdf("./CHP134.RARA.distribution.pdf", width=12, height=5)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

###########################################################################################CHP-134 SMAD4

files=Sys.glob("./data/1_CHP134/peaks/SMAD4/*.macs2.filter.narrowPeak")

files2=c()
for (f in files){
  #print(f)
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
files2

peaks <- lapply(files2,readPeakFile)
names(peaks)=c("DMSO_1","DMSO_2","RA_1","RA_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)

p2=plotDistToTSS(peakAnnoList,
                  title="Distribution of transcription factor-binding loci\nrelative to TSS")

pdf("./CHP134.SMAD4.distribution.pdf", width=12, height=5)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
dev.off()

#############################################################################################TGW-RARA
files=Sys.glob("./data/2_TGW/peaks/RARA/*.filter.narrowPeak")
files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
files2
peaks <- lapply(files2,readPeakFile)
names(peaks)=c("0_1","0_2","1_1","1_2","3_1","3_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)

p2=plotDistToTSS(peakAnnoList,
                  title="Distribution of transcription factor-binding loci\nrelative to TSS")

pdf("./TGW.RARA.distribution.pdf", width=12, height=5)

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

dev.off()

##########################################################################################TGW-SMAD4

files=Sys.glob("./data/2_TGW/peaks/SMAD4/*.filter.narrowPeak")

files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
files2
peaks <- lapply(files2,readPeakFile)
names(peaks)=c("0_1","0_2","1_1","1_2","3_1","3_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)

p2=plotDistToTSS(peakAnnoList,
                 title="Distribution of transcription factor-binding loci\nrelative to TSS")

pdf("./TGW.SMAD4.distribution.pdf", width=12, height=5)

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

dev.off()

##########################################################################################BE2-RARA
files=Sys.glob("./data/3_BE2/peaks/RARA/*.filter.narrowPeak")

files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
files2
peaks <- lapply(files2,readPeakFile)
names(peaks)=c("0_1","0_2","1_1","1_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)

p2=plotDistToTSS(peakAnnoList,
                 title="Distribution of transcription factor-binding loci\nrelative to TSS")

pdf("./BE2.RARA.distribution.pdf", width=12, height=5)

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

dev.off()
##########################################################################################BE2-SMAD4
files=Sys.glob("./data/3_BE2/peaks/SMAD4/*.filter.narrowPeak")

files2=c()
for (f in files){
  if(!str_detect(f,"FDR50")){
    files2=append(files2,f)
  }
}
files2
peaks <- lapply(files2,readPeakFile)
names(peaks)=c("0_1","0_2","1_1","1_2")
# peak

###################annotate peak with gene name and location information

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
###############genomic annotation
p1=plotAnnoBar(peakAnnoList)

p2=plotDistToTSS(peakAnnoList,
                 title="Distribution of transcription factor-binding loci\nrelative to TSS")

pdf("./BE2.SMAD4.distribution.pdf", width=12, height=5)

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

dev.off()