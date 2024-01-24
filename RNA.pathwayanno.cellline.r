setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/")
library("ggplot2")
library(tidyverse)
library("clusterProfiler")
library("enrichplot")
library("ggupset")
library("ggnewscale")
library(UpSetR)
library(VennDiagram)
library("ggridges")
library(RColorBrewer)
library(circlize)
library(gplots)
library(ggbiplot)
library(devtools)
library(ggbiplot)
library(genefilter)
library(plotly)
library('made4')
library(cowplot)
library("pasilla")
library("DESeq2")
library(dplyr)
library(lmtest)
library("ChIPseeker")
library("ChIPpeakAnno")
library(pheatmap)

#browseVignettes("ChIPpeakAnno")


################1. check the pathway enrichment for differential genes of each group and the overlapped core genes
mainDir <- "Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/enrichment_results_fpkm"
dir.create(file.path(mainDir))
########################input, read gmt files for Kegg,Go,Bioplanet and DSigDB, TFppi and histon factor
kegg=read.gmt("./list/enrichr/KEGG_2021_Human.txt")
go_b=read.gmt("./list/enrichr/GO_Biological_Process_2021.txt")
bplan=read.gmt("./list/enrichr/BioPlanet_2019.txt")
dsig=read.gmt("./list/enrichr/DSigDB.txt")
ppi=read.gmt("./list/enrichr/TRRUST_Transcription_Factors_2019.txt")
his=read.gmt("./list/enrichr/ENCODE_Histone_Modifications_2015.txt")
gmt=list(kegg=kegg,go=go_b,bplan=bplan,drug=dsig,tfppi=ppi,histone=his)

#######################annotated DEG genes for each group, to generating some basic results (plot&tables)
############Section 1 gene enrichment analysis
###########################################################################################
#browseVignettes("clusterProfiler")

files=Sys.glob("./data/rna/NBCellLines_2377244/*_FPKM.tsv")
for (f in files){
  s=read.table(f,header = T)
  deg=unique(s[s$P.Value<0.05&abs(s$log2FC)>1,2:3])
  all=unique(s[,2:3])
  all_or=all[order(-all$log2FC),]
  genelist=all_or[,2]
  names(genelist)=all_or[,1]
  deglist=deg[,1]
  l=length(deglist) ####DEG gene list
  a=length(genelist)####all genes in the table
  gname=str_extract(f,"(?<=STRANDED\\_).*(?=\\_FPKM)") ###extract group name
  subDir=gname
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
  dflist=data.frame(deglist)
  write.table(dflist, file = paste0(gname,"DEG.p0.05fc2.csv"), row.names=FALSE,col.names = F)##################output the DEG genelist
  ##################################generate and output enrichment analysis results for each group
  c=1
  for (db in gmt){
     x=enricher(deglist,TERM2GENE =db,pAdjustMethod="BH",qvalueCutoff = 0.25,pvalueCutoff = 0.5)
     
     y <- GSEA(genelist, TERM2GENE = db,pvalueCutoff = 0.05)
     dbname=names(gmt[c])
     print(dbname)
     dfx=data.frame(x)
     dfy=data.frame(y)
     if(length(dfx[,1])>=1){
       write.csv(dfx, file = paste0(dbname,".DEGenricher.csv"), row.names=FALSE)
       p1=barplot(x, showCategory=20)
       p2=dotplot(x, showCategory=20)
       p3=upsetplot(x)
       plotfile=paste0("./",dbname,".DEGenricher.pdf")
       pdf(plotfile, width=15, height=12)
       print(p1)
       print(p2)
       print(p3)
       dev.off()
       
     }
     if(length(dfy[,1])>=1){
       write.csv(dfy, file = paste0(dbname,".GSEA.csv"), row.names=FALSE)
       p2=dotplot(y, showCategory=20)
       p3=upsetplot(y)
       p4=ridgeplot(y)
      #p5=gseaplot2(y, geneSetID = 1:4, pvalue_table = TRUE)

       plotfile=paste0("./",dbname,".GSEA.pdf")
       pdf(plotfile, width=15, height=12)
       print(p2)
       print(p3)
       print(p4)
       dev.off()
     }
     
     c=c+1
     #break
  }
  #break
}

# pdf("./runningscore.GSEA.bplan.pdf", width=15, height=12)
# p5
# dev.off()

################################################################################################################
###############Section 2 PCA test, check the RNA seq sample quality
################################################################################################################
###################extract CPM value from files
setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/")
c=read.table("./data/rna/NBCellLines_2377244/GEELE-NBCellLines-2377244-STRANDED_RSEM_gene_count.2022-03-12_03-56-57.txt",header = T)
fk=read.table("./data/rna/NBCellLines_2377244/GEELE-NBCellLines-2377244-STRANDED_RSEM_gene_FPKM.2022-03-12_03-56-57.txt",header = T)
#########calculate colume sum, and get CMP
cmp=cbind(c[,1:4],log2(c[,c(5:28)]*1e6/colSums(c[,c(5:28)])))
fk=cbind(c[,1:4],log2(fk[,c(5:28)]))
#ex=distinct(cmp,geneSymbol,.keep_all = T)
ex=distinct(fk,geneSymbol,.keep_all = T)

sub2=ex[,c(5:28)]
row.names(sub2)=ex[,2]
sub2[sub2==0]=0.001
ins3=t(scale(t(sub2)))
ins4=ins3[complete.cases(ins3),]
length(ins4[,1])

#ins4$var=rowVars(ins4)
#s2=ins4[order(ins4$var,decreasing = T),]
#s3=subset (s2[1:3000,],select=-c(var))
pca=prcomp(t(ins4),scale=T,rank.=3)

tot <- summary(pca)[["importance"]]['Proportion of Variance',1:3]
tot_ratio <- 100 * sum(tot)

tit = paste0('Total Explained Variance =', tot_ratio )
pc1 = paste0('PC1 =', 100 * tot[1] )
pc2 = paste0('PC2 =', 100 * tot[2] )
pc3 = paste0('PC3 =', 100 * tot[3] )
names=str_extract(row.names(pca$x),"(?<=CHP.134.).*(?=\\.\\d)")
pcs <-data.frame(pca$x, group = names)
pcs$group=as.factor(pcs$group)


p1=ggplot(pcs,aes(x=PC1,y=PC2,color=group,label=group))+
   #geom_point(size=3)+
  geom_text(size=4)+
   ggtitle(tit)+
   xlab(pc1)+
   ylab(pc2)+
  theme_classic()

p2=ggplot(pcs,aes(x=PC1,y=PC3,color=group,label=group))+
  #geom_point(size=3)+
  geom_text(size=4)+
  #ggtitle(tit)+
  xlab(pc1)+
  ylab(pc3)+
  theme_classic()

p3=ggplot(pcs,aes(x=PC2,y=PC3,color=group,label=group))+
  #geom_point(size=3)+
  geom_text(size=4)+
  #ggtitle(tit)+
  xlab(pc2)+
  ylab(pc3)+
  theme_classic()
#######################heatmap


pdf("./total.pca.pdf", width=25, height=8)
plot_grid(p1,p2,p3,ncol=3)

dev.off()

pdf("./total.cluster.pdf", width=8, height=25)
heatplot(ins4)
dev.off()


