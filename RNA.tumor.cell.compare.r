setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/7_scripts")

library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library("clusterProfiler")
library("enrichplot")
library(cowplot)
####################################################Usage Description
##this script is used to compare gene expression level among  NB cells implanted to mouse (ctr vs drug) and CHP-134 cell line, and remove all other databases
##in v2.0, change the boxplot to stripchart
#############################################################################################################################################################
###############collecting RNA-seq data from different database
##NBcell line and xenograft mouse

d1=read.table("./data/5_RNA_tumor_cell_compare/CHP-134-tumor.STAR.RSEM.TPM.txt",header=T)
d2=read.table("./data/5_RNA_tumor_cell_compare/GEELE-NBCellLines-2377244-STRANDED_RSEM_gene_TPM.2022-03-12_03-56-57.txt",header = T)[,c(1,2,3,5,6,11,12,17,18,23,24)]
colnames(d2)=c("GeneID","geneSymbol","bioType","Ctl0d-1","Ctl0d-2","Ctl1d-1","Ctl1d-2","Ctl3d-1","Ctl3d-2","Ctl6d-1","Ctl6d-2")

d_cell_pdx=inner_join(d1,d2,by= c("GeneID","geneSymbol"))  %>%
  filter(bioType=='protein_coding')%>%
  group_by(geneSymbol)%>%
  filter(length(geneSymbol)==1)


########################################make the plot just for ID genes for the main figure
#####just use IDs and SMAD1
bmp_focus=read.table("./data/5_RNA_tumor_cell_compare/focus.bmp.signaling.kegg.list",header=T)

bmp_cell_pdx_focus = inner_join(bmp_focus,d_cell_pdx,by=c("GeneID"="geneSymbol"))%>%
  rowwise()%>%
  mutate(Cell_0=mean(c(`Ctl0d-1`,`Ctl0d-2`)),Cell_6=mean(c(`Ctl6d-1`,`Ctl6d-2`)))%>%
  mutate(PDX_base=median(c(`L49`,`L55`,`L59`,`L60`,`L73`)),
         PDX_drug=median(c(`L54`,`L56`,`L61`,`L72`)))%>%
  select(GeneID,Cell_0,Cell_6,PDX_base,PDX_drug)%>%
  pivot_longer(
    cols = Cell_0:PDX_drug,
    names_to = "type",
    values_to = "TPM"
  )

p1=ggplot(bmp_cell_pdx_focus,aes(x=type,y=GeneID,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 4,limits=c(0, 8))+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p1

###############id1 boxplot

sampleinfo=read.table("./data/5_RNA_tumor_cell_compare/sample.info-chp134.txt",header = T)
ID1_cell_pdx=t(d_cell_pdx[d_cell_pdx$geneSymbol=='ID1',-c(1,2,12,15,16,17,18)])%>%
  as.data.frame()

ID1_cell_pdx$sample=row.names(ID1_cell_pdx)
ID1_cell_pdx=inner_join(ID1_cell_pdx,sampleinfo,by="sample")
ID1_cell_pdx$data=paste(ID1_cell_pdx$source,ID1_cell_pdx$drug,sep="_")
ID1_cell_pdx=ID1_cell_pdx%>%
  select(V1,data)
colnames(ID1_cell_pdx)=c("ID1ex","data")


ID1_cell_pdx$data=factor(ID1_cell_pdx$data,levels=c("cellline_DMSO","cellline_RA","tumor_DMSO","tumor_RA"))

p2=ggplot(ID1_cell_pdx,aes(x=data,y=log2(ID1ex+1)))+
  geom_jitter(aes( color = data), 
              position = position_jitter(0.4), size = 8) +
  scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07","blue"))+
  stat_summary(color="black", size = 0.3,alpha=0.5,
               fun.data="mean_sdl",  fun.args = list(mult=1))+

  theme_classic(base_size = 20) +
  ggtitle(" ") +
  xlab("")+
  ylab("expression of ID1 (log2(TPM+1))")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")

p2

pdf("cellline.tumor.stripchart.heatmap.pdf", width=10, height=8)

plot_grid(p2, p1, nrow=1, rel_widths = c(1,1.5))

dev.off()






