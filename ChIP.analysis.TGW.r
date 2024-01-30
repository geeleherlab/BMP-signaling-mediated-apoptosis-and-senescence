
setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/7_scripts")

library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library(cowplot)
library(pheatmap)
library('made4')
library(PerformanceAnalytics)
library("clusterProfiler")
library("enrichplot")
library(corrplot)
library(R.utils)
library(ggvenn)
R.utils::setOption("clusterProfiler.download.method","auto")
#library("pathview")
#library('org.Hs.eg.db')
#####################catalog
##1. PCA/correlation
##2. differential chip and chip signal heatmap (commandline)
##3. GSEA function annotation
##4. Venn analysis to sum up overlaps between samples
##5. chip location annotation (gene, intergenic, exon,intron....)

##################################read files
saminfo=read.csv("./data/3_ChIP_GSEApathway_venn/2_TGW/sample.combined.chip.csv")
anno=read.table("./data/3_ChIP_GSEApathway_venn/2_TGW/RARA_SMAD4_combined.bed.anno",header = T)
d=read.table("./data/3_ChIP_GSEApathway_venn/2_TGW/all.counts.filter.data",header = T)
c_libsize=saminfo$libsize
###########################Get sample info
saminfo$source=factor(saminfo$source,levels = unique(saminfo$source))
tb=table(saminfo$source)
start=0
end=0
frag=list()
for(i in seq_along(tb)){
  start=end+1
  end=start+tb[i]-1
  frag[[i]]=c(start:end)
}
##########################################normalization peaks by library size and combined peaks intensity with peak location annotation
rownames(d)=d[,1]
d2=as.matrix(d[,-c(1)])
tshort=c_libsize/50000000
d3=t(t(d2)/tshort)
cutoff=30
keep <- rowSums(d3> cutoff) >= 2
d4 <- d3[keep, ]
df=data.frame(d4)
df$Region=rownames(df)
############################combining with gene annotation files
add_anno=inner_join(df,anno,by= "Region")[,-c(14:20)]

################## 1. PCA/correlation matrix 

df_pro=df[,-c(ncol(df))]
c=cor(df_pro)
pdf("./correlationTF.pdf", width=10, height=12)
testRes = cor.mtest(df_pro, conf.level = 0.95)
corrplot(c, p.mat = testRes$p,method = 'color',
         sig.level = c(0.001),
         type = 'lower',insig = 'blank', pch.col = 'grey20',
         addCoef.col ='black',diag=FALSE,
         col.lim = c(0.2, 1),
         col = COL2('PiYG')
)
dev.off()

#################################################pca
pca=prcomp(t(df_pro),scale=T,rank.=3)
tot <- summary(pca)[["importance"]]['Proportion of Variance',1:3]
tot_ratio <- 100 * sum(tot)

tit = paste0('Total Explained Variance =', tot_ratio )
pc1 = paste0('PC1 =', 100 * tot[1] )
pc2 = paste0('PC2 =', 100 * tot[2] )
pc3 = paste0('PC3 =', 100 * tot[3] )
names=str_extract(row.names(pca$x),"(?<=\\w\\.)\\w+$")
pcs <-data.frame(pca$x, group = names,label=row.names(pca$x))
pcs$group=as.factor(pcs$group)


p1=ggplot(pcs,aes(x=PC1,y=PC2,color=group,label=label))+
  #geom_point(size=3)+
  geom_text(size=4)+
  ggtitle(tit)+
  xlab(pc1)+
  ylab(pc2)+
  theme_classic()

p2=ggplot(pcs,aes(x=PC1,y=PC3,color=group,label=label))+
  #geom_point(size=3)+
  geom_text(size=4)+
  ggtitle(tit)+
  xlab(pc1)+
  ylab(pc3)+
  theme_classic()

pdf("./pca.TF.pdf", width=12, height=6)

plot_grid(p1, p2, label_size = 12)
dev.off()

##############################################################2. deferential peaks calling for each chip type. gene enrichment pathway annotation. GSEA annotation for deferential peaks

########################do differential peaks calling and adding to annotation table for each chip
for (i in seq_along(tb)){

    d3=as.matrix(add_anno[,frag[[i]][1:4]])
    mode(d3) <- "integer"
    #samra=saminfo[c(frag[[i]][1:2],tail(frag[[i]],n=2)),]
    samra=saminfo[frag[[i]][1:4],]
    samra$type="RA"
    samra[1:2,]$type="DMSO"
    samra$type=as.factor(samra$type)

    dds <- DESeqDataSetFromMatrix(countData = d3,
         colData = samra, design = ~ type)
   
    sizeFactors(dds)=c(1,1,1,1)
    ddsTC <- DESeq(dds)
    ddsTC$sizeFactor
    resultsNames(ddsTC)
    resTC <- results(ddsTC,name = resultsNames(ddsTC)[2],test='Wald')
    diff=data.frame(resTC@listData)
    colnames(diff)=str_c(unique(samra$source),colnames(diff),sep = "_")
    add_anno=cbind(add_anno,diff)
}


#######################differential peaks
add_diff=add_anno
for (i in seq_along(levels(samra$source))){
   name=paste(levels(samra$source)[i],"diffpeak",sep="_")
   fc=paste(levels(samra$source)[i],"log2FoldChange",sep="_")
   pv=paste(levels(samra$source)[i],"pvalue",sep="_")
   base=paste(levels(samra$source)[i],"baseMean",sep="_")
   plot=paste(levels(samra$source)[i],"plot",sep="_")
   add_diff[,name]="nochange"
   add_diff[add_diff[[fc]]>0 & !is.na(add_diff[[pv]]) & add_diff[[pv]]<0.05,name]="up"
   add_diff[add_diff[[fc]]<0 & !is.na(add_diff[[pv]]) & add_diff[[pv]]<0.05,name]="down"
   p=ggplot(data=add_diff[add_diff[[base]]>20,], aes(x=!!sym(fc), y=-log10(!!sym(pv)),color=!!sym(name))) +
     geom_point(alpha=0.3) +
     scale_color_manual(values=c("blue", "black", "red")) +
     theme_classic(base_size =20)+
     geom_hline(yintercept=-log10(0.05), col="red")
     filename=paste("./",name,".3.pdf",sep="")
     pdf(filename, width=8, height=12)
     print(p)
     dev.off()
}
#####################3 Annotation by GSEA and enricher and genelist
bplan=read.gmt("./data/3_ChIP_GSEApathway_venn/BioPlanet_2019.txt")
for (i in seq_along(levels(saminfo$source))){
  base=paste(levels(saminfo$source)[i],"baseMean",sep="_")
  diff=paste(levels(saminfo$source)[i],"diffpeak",sep="_")
  df=add_diff[order(add_diff[[base]],decreasing = T),]%>%
    filter(!!sym(base)>=mean(add_diff[[base]])*0.5)
  or_list=unique(df[,c("Closest_Gene",base)])
  or_list=na.omit(or_list)
  genelist=or_list[[base]]
  names(genelist)=or_list$Closest_Gene
  y <- GSEA(genelist, TERM2GENE = bplan,pvalueCutoff = 0.5)
  p1=gseaplot2(y, geneSetID = 1:4, pvalue_table = TRUE,base_size = 25)
  dfy=data.frame(y)
  d2=dfy %>%
    slice_min(pvalue,n=20)
  
  d2$Description=factor(d2$Description,levels = d2[order(d2$pvalue,decreasing = T),'Description'])
  p2=ggplot(d2,aes(x=Description,y=enrichmentScore,fill=-log10(pvalue)))+
    geom_bar(stat="identity")+
    scale_fill_continuous(low="blue", high="red")+
    coord_flip()+
    theme_classic(base_size = 25)
  
  filename=paste("./",levels(saminfo$source)[i],".gsea.all.csv",sep="")
  write.csv(dfy,filename,quote = F)
  filename=paste("./",levels(saminfo$source)[i],".gene.annotation.pdf",sep="")
  ####gene annotation plot
  pdf(filename, width=18, height=12)
  if(dim(dfy)>1){
    print(p1)
    print(p2)
  }
  
  degtype=c("up","down")
  for (j in seq_along(degtype) ){
    deglist=df%>%
      filter(!!sym(diff)==degtype[j]) %>%
      pull(Closest_Gene) %>%
      unique()
    dfx=data.frame()
    if(length(deglist)>30){
      x=enricher(deglist,TERM2GENE =bplan,pAdjustMethod="BH",qvalueCutoff = 0.05)
      dfx=data.frame(x)
      filename=paste("./",levels(saminfo$source)[i],".",degtype[j],".enricher.deseq.csv",sep="")
      write.csv(dfx,filename,quote = F)
      p3=dotplot(x, showCategory=16,font.size=25,title=filename)
      if(dim(dfx)>1){
        print(p3)
      }
    }
  }
  dev.off()
}

#################################################################4. Venn analysis to sum up overlaps between samples

peaks=paste("peaks",1:nrow(add_anno))
add_anno$peaks=peaks

x=list()
for (i in seq_along(tb)){
  base=paste(levels(saminfo$source)[i],"baseMean",sep="_")
  keep=add_anno[[base]]>median(d4)
  x[[levels(saminfo$source)[i]]]=add_anno[keep,]$peaks
}

p=ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 8,text_size=6
)

pdf("./venn.peaksoverlap.pdf", width=18, height=12)
p
dev.off()