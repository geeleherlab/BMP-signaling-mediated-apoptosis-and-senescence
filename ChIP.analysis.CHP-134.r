setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/7_scripts")

library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library(cowplot)
library("clusterProfiler")
library("enrichplot")
library(ggvenn)
#####################catalog
#1. differential chip calling
#2. gene annotation for the overlapped loci between RARA and Smad9 TFs 

##################################read files

saminfo=read.csv("./data/3_ChIP_GSEApathway_venn/1_CHP134/sample.combined.chip.csv")
anno=read.table("./data/3_ChIP_GSEApathway_venn/1_CHP134/combined.all.reprodu.bed.anno",header = T)
d=read.table("./data/3_ChIP_GSEApathway_venn/1_CHP134/all.counts.filter.data",header = T)
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
add_anno=inner_join(df,anno,by= "Region")[,-c(18:24)]

##############################################################1. deferential peaks calling for each chip type. gene enrichment pathway annotation. GSEA annotation for deferential peaks
 ########################do differential peaks calling and adding to annotation table for each chip
for (i in seq_along(tb)){

    #d3=as.matrix(add_anno[,c(frag[[i]][1:2],tail(frag[[i]],n=2))])
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
}

peaks=paste("peaks",1:nrow(add_anno))
add_diff$peaks=peaks

x=list()
for (i in 1:2){
  base=paste(levels(saminfo$source)[i],"baseMean",sep="_")
  keep=add_anno[[base]]> median(d4)
  x[[levels(saminfo$source)[i]]]=add_diff[keep,]$peaks
}

p=ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 8,text_size=6
)

pdf("./venn.peaks.RARA.SMAD9.overlap.pdf", width=18, height=12)
p
dev.off()

########################################2 Annotation for overlapped loci

a=Reduce(intersect,x)%>%
  data.frame(peaks=.)

over=inner_join(add_diff,a,by= "peaks")

#################################a baseline overlap annotation
bplan=read.gmt("./data/3_ChIP_GSEApathway_venn/BioPlanet_2019.txt")
base="SMAD9_baseMean"
df=over[order(over[[base]],decreasing = T), ]

or_list=unique(df[,c("Closest_Gene",base)])
or_list=na.omit(or_list)
genelist=or_list[[base]]
names(genelist)=or_list$Closest_Gene

y <- GSEA(genelist, TERM2GENE = bplan,pvalueCutoff = 0.2)
dfy1=data.frame(y)
View(dfy1)

p1=gseaplot2(y, geneSetID = 1:3, pvalue_table = TRUE,base_size = 25)

x=enricher(over$Closest_Gene,TERM2GENE =bplan,pAdjustMethod="BH",qvalueCutoff = 0.05)
dfx1=data.frame(x)

###################################b increased RARA + SMAD overlap

over_up=over[over$RARA_diffpeak == 'up',]

base="SMAD9_baseMean"
df=over_up[order(over_up[[base]],decreasing = T),]
or_list=unique(df[,c("Closest_Gene",base)])
or_list=na.omit(or_list)
genelist=or_list[[base]]
names(genelist)=or_list$Closest_Gene

y <- GSEA(genelist, TERM2GENE = bplan,pvalueCutoff = 0.2)
dfy=data.frame(y)
View(dfy)

p1=gseaplot2(y, geneSetID = 1:4, pvalue_table = TRUE,base_size = 25)

x=enricher(over_up$Closest_Gene,TERM2GENE =bplan,pAdjustMethod="BH",qvalueCutoff = 0.05)
dfx=data.frame(x)
View(dfx)
p3=dotplot(x, showCategory=16,font.size=25,title="overlop loci gene annotation")
d2=dfy %>%
  slice_min(pvalue,n=20)

d2$Description=factor(d2$Description,levels = d2[order(d2$pvalue,decreasing = T),'Description'])
p2=ggplot(d2,aes(x=Description,y=enrichmentScore,fill=-log10(pvalue)))+
  geom_bar(stat="identity")+
  scale_fill_continuous(low="blue", high="red")+
  coord_flip()+
  theme_classic(base_size = 25)

pdf("upRARA_RARASMAD9overlap.gene.annotation.pdf", width=18, height=12)
print(p1)
print(p2)
print(p3)
dev.off()

write.csv(dfx,"upRARA_RARASMAD9overlap.nnricer.anno.csv",quote = F)
write.csv(dfy,"upRARA_RARASMAD9overlap.gsea.anno.csv",quote = F)

####################################increased signal at pathway level
##box-plot to show chip signal at pathway level

##apoptosis, tgf-beta, averaged chip signal for all overlaped loci
apop_l=c()###
tgf_l=c()###

for (i in seq_along(dfy1$Description)){
  a=str_split(dfy1$core_enrichment[i],"\\/",simplify = T)%>%
    as.vector()
  if(str_detect(dfy1$Description[i],"TGF|BMP|ALK1")){

    tgf_l=c(tgf_l,a)
    
  }
  
  else if(str_detect(dfy1$Description[i],"apoptosis")){
    apop_l=c(apop_l,a)
    print(dfy1$Description[i])
  }
}

apop_df=unique(apop_l) %>%
  data.frame()%>%
  `colnames<-`(c("geneID")) %>%
  inner_join(over,.,by=c("Closest_Gene"="geneID"))%>%
  mutate(pathway="apoptosis") 

tgf_df=unique(tgf_l) %>%
  data.frame()%>%
  `colnames<-`(c("geneID")) %>%
  inner_join(over,.,by=c("Closest_Gene"="geneID"))%>%
  mutate(pathway="TGF-beta")
over$pathway="all"

#####################combine all data-set into one data-frame for plot baseline level comparison
all_df=rbind(apop_df,tgf_df,over)%>%
  pivot_longer(
    cols = c(RARA_baseMean,SMAD9_baseMean,SMAD4_baseMean),
    names_to = "Type",
    names_pattern = "(.*)_baseMean",
    values_to = "intensity"
  )%>%
  select(pathway,Type,intensity)

###Long format
p=ggplot(all_df,aes(x=pathway,y=intensity))+
  facet_grid(. ~ Type)+
  geom_violin(aes(color=pathway),width=0.7,lwd=1,scale="width")+
  geom_boxplot(width=0.1,aes(fill=pathway),outlier.shape = NA,lwd=0.8) + 
  theme_classic(base_size = 25)+
  ylim(0,300)+
  theme(legend.position="none")+
  ylab("ChIP binding intensity")+
  xlab("")

 pdf("chip-insensity-pathway.baseline.pdf", width=16, height=8)
 print(p) 
 dev.off()

############################showing changes between dmso and 1day treatment

all_df_tr=rbind(apop_df,tgf_df,over)%>%
  rowwise() %>%
  mutate(RARA_0_mean=mean(c(DMSO_1.RARA, DMSO_2.RARA)), RARA_1_mean=mean(c(RA_1.RARA, RA_2.RARA)),
        SMAD9_0_mean=mean(c(Min_0_1.SMAD9, Min_0_2.SMAD9)), SMAD9_1_mean=mean(c(Min_1_1.SMAD9, Min_1D_2.SMAD9)),
        SMAD9_3_mean=mean(c(Min_Day3_1.SMAD9, Min_Day3_2.SMAD9)), SMAD9_6_mean=mean(c(Min_6_1.SMAD9, Min_6_2.SMAD9)),
        SMAD4_0_mean=mean(c(DMSO.Rep1.SMAD4, DMSO.Rep2.SMAD4)), SMAD4_1_mean=mean(c(ATRA.Rep1.SMAD4, ATRA.Rep2.SMAD4)) 
        )%>%
  pivot_longer(
    cols = c(RARA_0_mean:SMAD4_1_mean),
    names_to = c("Type","Time"),
    names_pattern = "(\\w+)_(\\d+)_mean",
    values_to = "intensity"
  )%>%
  select(pathway,Type,Time,intensity)
 
p=ggplot(all_df_tr,aes(x=pathway,y=intensity,fill=Time))+
  facet_grid(. ~ Type)+
 # geom_violin(width=0.5,lwd=1,scale="width")+
  geom_boxplot(width=0.6,outlier.shape = NA,lwd=0.8) + 
  theme_classic(base_size = 25)+
  ylim(0,300)+
  #theme(legend.position="none")+
  ylab("ChIP binding intensity")+
  xlab("")

pdf("chip-insensity-pathway.add.times.pdf", width=16, height=8)
print(p) 
dev.off()
 