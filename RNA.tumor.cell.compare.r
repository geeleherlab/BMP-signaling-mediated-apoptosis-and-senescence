setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/3_RNA/")
library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library("clusterProfiler")
library("enrichplot")
library(cowplot)
##########################################workspace
saminfo=read.table("./data/sample.info-chp134.txt",header = T)
#samcontrol=saminfo[saminfo$drug=='DMSO',]
#samcontrol$source=as.factor(samcontrol$source)

##################get TPM
d1=read.table("./data/CHP-134-tumor.STAR.RSEM.TPM.txt",header=T)
d2=read.table("./data/GEELE-NBCellLines-2377244-STRANDED_RSEM_gene_TPM.2022-03-12_03-56-57.txt",header = T)[,c(1,2,3,5,6,11,12,17,18,23,24)]
colnames(d2)=c("GeneID","geneSymbol","bioType","Ctrl-0d-1","Ctrl-0d-2","Ctrl-1d-1","Ctrl-1d-2","Ctrl-3d-1","Ctrl-3d-2","Ctrl-6d-1","Ctrl-6d-2")
c_gexp=inner_join(d1,d2,by= c("GeneID","geneSymbol"))

#######LOAD TCGA and NB patient tumor expression data

load("../../1_SMAD9_survival/dataIn/allExprData.RData")
load("../../1_SMAD9_survival/dataIn/tenRuvNewStandardApproach.RData")
nb <- read.delim("../../1_SMAD9_survival/list/PanNBL_Sig18_Expression.diagonly.trusSeqOnly.utest.log.CtoA.show.noFilter.mycn.17q.txt", as.is=T)


dt=as.data.frame(c_gexp %>%
                   filter(bioType=='protein_coding') %>%
                   group_by(geneSymbol)%>%
                   filter(length(geneSymbol)==1)%>%
                   select(-bioType))

target=c("ID1","ID2","ID3","ID4")

for (i in seq_along(target)){
  
  d_s=as.data.frame(t(dt[dt$geneSymbol==target[i],c(-1,-2)]))
  d_s$sample=row.names(d_s)
  colnames(d_s)=c("TPM","sample")
  
  dtpm=inner_join(d_s,saminfo,by="sample")%>%
    filter(drug=='DMSO')
  
  dtpm$source=as.factor(dtpm$source)
  
  tcga_expr=data.frame(ex=log2(exp(allExprData[target[i],])),type="TCGA")
  nb_expr=data.frame(ex=log2(as.numeric(nb[nb$GeneSet==target[i],-c(1,90,91)])+1),type="NB-patient")
  
  #################TCGA, NBpatient, CHP-134 cells implant in mouse and CHP-134 cell in vitro expression comparison
  
  mtm_expr=data.frame(ex=log2(dtpm[dtpm$source=='tumor' & dtpm$drug=='DMSO','TPM']),type="cell-in-Mouse")
  cell_expr=data.frame(ex=log2(dtpm[dtpm$source=='cellline' & dtpm$drug=='DMSO','TPM']),type="cell-in-vitro")
  comb=rbind(tcga_expr,nb_expr,mtm_expr,cell_expr)
  comb$type=as.factor(comb$type)
  
  p1=ggplot(data=comb, aes(x=type, y=ex,fill=type)) +
    geom_boxplot(alpha=0.6)+
    ylab("log2(TPM)")+
    theme_classic(base_size = 30)+
    xlab("")+
    theme(axis.text.x=element_text(angle=90, hjust=1))

  
  pdf(paste(target[i],".expression.pdf",sep=""), width=10, height=8)
  
  print(p1)
  
  dev.off()
  
  print(target[i])
  
}

################combined tumor and cell line gene expression table
d1=read.table("./data/CHP-134-tumor.STAR.RSEM.gene.count.txt",header=T)
d2=read.table("./data/GEELE-NBCellLines-2377244-STRANDED_RSEM_gene_count.2022-03-12_03-56-57.txt",header = T)[,c(1,2,3,5,6,11,12,17,18,23,24)]
colnames(d2)=c("GeneID","geneSymbol","bioType","Ctrl-0d-1","Ctrl-0d-2","Ctrl-1d-1","Ctrl-1d-2","Ctrl-3d-1","Ctrl-3d-2","Ctrl-6d-1","Ctrl-6d-2")
c_gexp=inner_join(d1,d2,by= c("GeneID","geneSymbol"))

##################filter genes with duplicate name and not belongs to protein coding genes.....

d3=as.data.frame(c_gexp %>%
                   filter(bioType=='protein_coding') %>%
                   group_by(geneSymbol)%>%
                   filter(length(geneSymbol)==1))


row.names(d3)=d3[,2]
d4=as.matrix(d3[,!names(d3) %in% c("GeneID","geneSymbol","bioType")])
mode(d4) <- "integer"

######################################create samples info and count matrix for each comparison
###considering 5 different comparison groups
c1=which(saminfo$source == "tumor")
c2=which(saminfo$drug=="DMSO")
c3=which((saminfo$source=="tumor"&saminfo$drug=="RA")|(saminfo$source=="cellline"&saminfo$drug=="DMSO"))
c4=which(saminfo$source == "cellline"&(saminfo$time==0|saminfo$time==1))
c5=which(saminfo$source == "cellline"&(saminfo$time==0|saminfo$time==3))
c6=which(saminfo$source == "cellline"&(saminfo$time==0|saminfo$time==6))

tb=list(c1,c2,c3,c4,c5,c6)
names(tb)=c("Tumor_Bvsdrug","BcellvsBTumor","BcellvsDrugTumor","Cell_Bvs1d","Cell_Bvs3d","Cell_Bvs6d")

add_diff=cbind(d3[,1:2],d4)

for (i in seq_along(tb)){
    ds=d4[,tb[[i]]]
    samra=saminfo[tb[[i]],]
    if(i==2|i==3){

      samra$source=as.factor(samra$source)
      dds <- DESeqDataSetFromMatrix(countData = ds,
                                    colData = samra,
                                    design = ~ source)
      
    }
    
    else {
      samra$drug=as.factor(samra$drug)
      dds <- DESeqDataSetFromMatrix(countData = ds,
                                    colData = samra,
                                    design = ~ drug)
    }

    ddsTC <- DESeq(dds)
    resultsNames(ddsTC)
    resTC <- results(ddsTC,name = resultsNames(ddsTC)[2],test='Wald',independentFiltering=FALSE,cooksCutoff=FALSE)
    diff=data.frame(resTC@listData)
    
    diff$diffpeak="nochange"
    diff[is.na(diff$log2FoldChange)==F& diff$log2FoldChange>1 & is.na(diff$padj)==F & diff$padj<0.05,"diffpeak" ]="up"

    diff[is.na(diff$log2FoldChange)==F& diff$log2FoldChange< -1 & is.na(diff$padj)==F & diff$padj<0.05,"diffpeak" ]="down"
    
    
    p=ggplot(data=diff, aes(x=log2FoldChange, y=-log10(padj),color=diffpeak)) +
      geom_point(alpha=0.3) +
      scale_color_manual(values=c("blue", "black", "red")) +
      theme_classic(base_size =20)
    #geom_hline(yintercept=-log10(0.05), col="red")
    p
    filename=paste("./",names(tb)[i],".differential.genes.pdf")
    
    pdf(filename, width=8, height=8)
    print(p)
    dev.off()
    
    colnames(diff)=str_c(names(tb)[i],colnames(diff),sep = "_")
    add_diff=cbind(add_diff,diff)
    
}

write.csv(add_diff,"./all.rawcounts.Deseq2.metrics.csv",quote=F)

table(add_diff$Tumor_Bvsdrug_diffpeak)
#############################################################################################correlation analysis
library(ggpubr)
library(ggcyto)
library(flowViz)

col1 <- densCols(add_diff$BcellvsBTumor_log2FoldChange, add_diff$Cell_Bvs1d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))

p1=ggplot(data=add_diff, aes(y=BcellvsBTumor_log2FoldChange, x=Cell_Bvs1d_log2FoldChange)) +
  geom_point(color=col1,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p1

col2 <- densCols(add_diff$BcellvsBTumor_log2FoldChange, add_diff$Cell_Bvs3d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))

p2=ggplot(data=add_diff, aes(y=BcellvsBTumor_log2FoldChange, x=Cell_Bvs3d_log2FoldChange)) +
  geom_point(color=col2,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p2


col3 <- densCols(add_diff$BcellvsBTumor_log2FoldChange, add_diff$Cell_Bvs6d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))
p3=ggplot(data=add_diff, aes(y=BcellvsBTumor_log2FoldChange, x=Cell_Bvs6d_log2FoldChange)) +
  geom_point(color=col3,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p3

col4 <- densCols(add_diff$BcellvsDrugTumor_log2FoldChange, add_diff$Cell_Bvs1d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))
p4=ggplot(data=add_diff, aes(y=BcellvsDrugTumor_log2FoldChange, x=Cell_Bvs1d_log2FoldChange)) +
  geom_point(color=col4,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p4

col5 <- densCols(add_diff$BcellvsDrugTumor_log2FoldChange, add_diff$Cell_Bvs3d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))
p5=ggplot(data=add_diff, aes(y=BcellvsDrugTumor_log2FoldChange, x=Cell_Bvs3d_log2FoldChange)) +
  geom_point(color=col5,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p5

col6 <- densCols(add_diff$BcellvsDrugTumor_log2FoldChange, add_diff$Cell_Bvs6d_log2FoldChange, colramp = flowViz::flowViz.par.get("argcolramp"))
p6=ggplot(data=add_diff, aes(y=BcellvsDrugTumor_log2FoldChange, x=Cell_Bvs6d_log2FoldChange)) +
  geom_point(color=col6,size=0.5) +
  geom_smooth(method = "lm", se = FALSE)+
  theme_classic(base_size = 20)+
  stat_cor(method = "pearson",p.accuracy = 1e-3, r.accuracy = 0.01,label.y=15,label.x=-5)
p6


pdf("./tumorbase_celline_stages.pdf", width=25, height=18)
plot_grid(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()

pdf("./tumorbase_celline_stages.justbase.pdf", width=25, height=9)
plot_grid(p1,p2,p3,nrow = 1)
dev.off()

#############################################################################################GSEA annotation

bplan=read.gmt("../list/enrichr/BioPlanet_2019.txt")

order_fc=add_diff[order(add_diff$BcellvsBTumor_log2FoldChange,decreasing = T),] %>%
  dplyr::filter(!is.na(BcellvsBTumor_log2FoldChange))

genelist=order_fc$BcellvsBTumor_log2FoldChange
names(genelist)=order_fc$geneSymbol

y <- GSEA(genelist, TERM2GENE = bplan,pvalueCutoff = 0.4)
p=gseaplot2(y, geneSetID = 6:10, pvalue_table = TRUE,base_size = 25)
dfy=data.frame(y)
View(dfy)
write.csv(dfy,"./BcellvsBTumor_log2FoldChange.GSEA.annotation.csv",quote=F)

p3=gseaplot2(y, geneSetID = c(11,43,68,97,122), pvalue_table = TRUE,base_size = 25)

d2=data.frame()
for (i in seq_along(dfy$Description)){
  if(str_detect(dfy$Description[i],"TGF|BMP|SMAD|Inhibitor|Retinol|Cell cycle|cell cycle|G0|G1|G2|Axon guidance|Tyrosine|DNA replication")& dfy$pvalue<0.05){
    #if(dfy$pvalue[i] < 0.05){
      d2=rbind(d2,dfy[i,])
    #}
  }
}


d2$Description=factor(d2$Description,levels = d2[order(d2$pvalue,decreasing = T),'Description'])

p1=ggplot(d2,aes(x=Description,y=enrichmentScore,fill=-log10(pvalue)))+
  geom_bar(stat="identity")+
  scale_fill_continuous(low="blue", high="red")+
  coord_flip()+
  theme_classic(base_size = 20)+
  xlab("enriched pathway for tumor from implanted cell")

#####################cell line annotation
order_fc=add_diff[order(add_diff$Cell_Bvs6d_log2FoldChange,decreasing = T),] %>%
  dplyr::filter(!is.na(Cell_Bvs6d_log2FoldChange))

genelist=order_fc$Cell_Bvs6d_log2FoldChange
names(genelist)=order_fc$geneSymbol

y2 <- GSEA(genelist, TERM2GENE = bplan,pvalueCutoff = 0.4)

dfy2=data.frame(y2)
View(dfy2)

p4=gseaplot2(y2, geneSetID = c(10,45,140,226,76), pvalue_table = TRUE,base_size = 25)

write.csv(dfy2,"./Bcellvs6daycell_log2FoldChange.GSEA.annotation.csv",quote=F)

#d_cell=read.csv("./chp134-tumor/cellline6dvs0d.kegg.GSEA.csv",header = T)
dc=data.frame()
for (i in seq_along(dfy2$Description)){
  if(str_detect(dfy2$Description[i],"TGF|BMP|SMAD|Inhibitor|Retinol|Cell cycle|cell cycle|G0|G1|G2|Axon guidance|Tyrosine|DNA replication")){
    #if(dfy2$pvalue[i] < 0.05){
      dc=rbind(dc,dfy2[i,]) 
    #}
  }
}

dc$Description=factor(dc$Description,levels = dc[order(dc$pvalue,decreasing = T),'Description'])


p2=ggplot(dc,aes(x=Description,y=enrichmentScore,fill=-log10(pvalue)))+
  geom_bar(stat="identity")+
  scale_fill_continuous(low="blue", high="red")+
  coord_flip()+
  theme_classic(base_size = 20)+
  xlab("enriched pathway for 6 days RA-treated cell in vitro")

pdf("./BcellvsBTumor_Bcellvscell6day.select.GSEA.barplot.pdf", width=15, height=20)
plot_grid(p1,p2,nrow = 2,rel_heights  = c(2, 2))
dev.off()


pdf("./BcellvsBTumor_Bcellvscell6day.select.GSEA.lineplot.pdf", width=24, height=12)
plot_grid(p3,p4,nrow = 1,rel_heights  = c(2, 2))
dev.off()
