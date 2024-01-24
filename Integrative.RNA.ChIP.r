setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/4_integrate/")
#setwd( "/Volumes/groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/chip_consensus_ref/" )
library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library(cowplot)
library(ggridges)
##################################content
##this scripts is trying to integrate information from both RNA-seq and ChIP-seq data.
##1. first, calculate the overlaps between differential expressed genes and differential changed RARA peaks 
##2. then, for genes that involved in apoptosis, cell cycle, differentiation and TGF signalling pathway, 
##exploring how much of them are DEG and their changes in 1,3 and 6 RA treatment. and the distribution of their foldchanges

##################################read files
###read RNA and ChIP information
###read chip peaks overlapped between RARA and SMAD4, which also include information of differential peaks information for two TFs
dc=read.csv("../1_ChIP/results/1_CHP134/RARA_SMAD4/overlap.SMAD4.RARA.csv")[,c("Closest_Gene","Distance","RARA_diffpeak","SMAD4_diffpeak")] %>%
  group_by(Closest_Gene) %>%
  filter(Distance==min(Distance)) %>%
  distinct(Closest_Gene, .keep_all = TRUE) %>%
  ungroup()

###RNA-seq
dr=read.csv("../3_RNA/all.rawcounts.Deseq2.metrics.csv")

#################################Main:
###1. find the overlaps between differential expressed genes and differential changed RARA peaks
cb=left_join(dr,dc,by= c("geneSymbol"="Closest_Gene"))%>% 
  replace_na(list(RARA_diffpeak = 'no_peak', SMAD4_diffpeak = 'no_peak')) %>%
  select(Cell_Bvs1d_diffpeak,Cell_Bvs3d_diffpeak,Cell_Bvs6d_diffpeak,RARA_diffpeak,SMAD4_diffpeak)

cb_l<-cb %>%
  pivot_longer(
    cols = c(Cell_Bvs1d_diffpeak,Cell_Bvs3d_diffpeak,Cell_Bvs6d_diffpeak),
    names_to = c("Time"),
    names_pattern = "(.*)_diffpeak",
    values_to = "diffexpr"
  )
### sum up the counts for the overlaps between differential genes and differential chip peaks
sum_rara=cb_l %>%
  group_by(Time,diffexpr,RARA_diffpeak) %>%
  summarise(n=n())

sum_smad4=cb_l %>%
  group_by(Time,diffexpr,SMAD4_diffpeak) %>%
  summarise(n=n())

table(cb$Cell_Bvs3d_diffpeak,cb$SMAD4_diffpeak)

#############plot
sum_rara$RARA_diffpeak=factor(sum_rara$RARA_diffpeak,levels=c("no_peak","nochange","up","down"))
sum_smad4$SMAD4_diffpeak=factor(sum_smad4$SMAD4_diffpeak,levels=c("no_peak","nochange","up","down"))

colorvalue=c("grey","#FFDE7D","#F6416C","#00B8A9")
p1=ggplot(sum_rara[sum_rara$diffexpr!='nochange',],aes(x=diffexpr,y=n,fill=RARA_diffpeak))+
  scale_fill_manual(values = colorvalue)+
  facet_grid(. ~ Time)+
  geom_bar(stat="identity")+
  theme_classic(base_size = 25)+
  #theme(legend.position="none")+
  ylab("the number of DEG")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))
p1

p2=ggplot(sum_rara[sum_rara$diffexpr=='nochange',],aes(x=diffexpr,y=n,fill=RARA_diffpeak))+
  scale_fill_manual(values = colorvalue)+
  facet_grid(. ~ Time)+
  geom_bar(stat="identity")+
  theme_classic(base_size = 25)+
  #theme(legend.position="none")+
  ylab("the number of DEG")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))
p2

p3=ggplot(sum_smad4[sum_smad4$diffexpr!='nochange',],aes(x=diffexpr,y=n,fill=SMAD4_diffpeak))+
  scale_fill_manual(values = colorvalue)+
  facet_grid(. ~ Time)+
  geom_bar(stat="identity")+
  theme_classic(base_size = 25)+
  #theme(legend.position="none")+
  ylab("the number of DEG")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))
p3

p4=ggplot(sum_smad4[sum_smad4$diffexpr=='nochange',],aes(x=diffexpr,y=n,fill=SMAD4_diffpeak))+
  scale_fill_manual(values = colorvalue)+
  facet_grid(. ~ Time)+
  geom_bar(stat="identity")+
  theme_classic(base_size = 25)+
  #theme(legend.position="none")+
  ylab("the number of DEG")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))
p4

##################output
pdf("overlap.DEG.ChIP.statistic2.pdf", width=19, height=12)
plot_grid(p1, p2, p3,p4, nrow = 2, label_size = 12)
dev.off()




###2. how much of apoptosis, cell cycle, differentiation and TGF involved genes are DEG and their changes in 1,3 and 6 RA treatment. and the distribution of their foldchanges

dp=read.csv("./data/genes.in.four.pathway.addchip.csv")%>%
  filter(cellline =='CHP134')####genes in pathway, and also have chip signal of RARA&SMAD4
dp$pathway=as.factor(dp$pathway)

df=data.frame()

for (i in seq_along(levels(dp$pathway))){
  a=left_join(dp[dp$pathway==levels(dp$pathway)[i],],dr,by=c("Closest_Gene"="geneSymbol")) %>%
    select(Closest_Gene,pathway,Cell_Bvs1d_diffpeak,Cell_Bvs1d_log2FoldChange,Cell_Bvs3d_diffpeak,Cell_Bvs3d_log2FoldChange,Cell_Bvs6d_diffpeak,Cell_Bvs6d_log2FoldChange )
  df=rbind(df,a)
}
head(df)

df_l<-df %>%
  pivot_longer(
    cols = c(Cell_Bvs1d_diffpeak,Cell_Bvs3d_diffpeak,Cell_Bvs6d_diffpeak),
    names_to = c("Time"),
    names_pattern = "(.*)_diffpeak",
    values_to = "diffexpr"
  )

df_v<-df %>%
  pivot_longer(
    cols = c(Cell_Bvs1d_log2FoldChange,Cell_Bvs3d_log2FoldChange,Cell_Bvs6d_log2FoldChange),
    names_to = c("Time"),
    names_pattern = "(.*)_log2FoldChange",
    values_to = "foldchange"
  )

sum_path=df_l %>%
  group_by(Time,pathway,diffexpr) %>%
  summarise(n=n())


sum_path$diffexpr=factor(sum_path$diffexpr,levels=c("nochange","up","down"))
p1=ggplot(sum_path[sum_path$diffexpr!="nochange",],aes(x=pathway,y=n,fill=diffexpr))+
  scale_fill_manual(values = c("#F6416C","#19A7CE"))+
  facet_grid(. ~ Time)+
  geom_bar(stat="identity")+
  theme_classic(base_size = 20)+
  #theme(legend.position="none")+
  ylab("the number of DEG")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))
p1

p2=ggplot(df_v, aes(x=pathway,y=foldchange,fill=pathway))+
  #scale_fill_manual(values = c("#F6416C","#00B8A9"))+
  facet_grid(. ~ Time)+
  geom_boxplot()+
  theme_classic(base_size = 20)+
  coord_cartesian( ylim = c(-2, 2))+
  ylab("log2FoldChange distribution")+
  xlab("")+
  theme(strip.background.x=element_rect(color = NA,  fill=NA))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "#F7C04A", size=1.5)
p2

pdf("pathway,expression.pdf", width=18, height=6)
plot_grid(p1, p2,  nrow = 1, label_size = 12)
dev.off()


