
setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/7_scripts")
library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library("clusterProfiler")
library("enrichplot")
library(cowplot)
###############################################################################################Usage Description
##this script is used to read interested gene expression from the programs annotations from single cell rnaseq data of 5 database
##and combined them together into a boxplot
###############################################################################################
t <- Sys.glob("./data/8_Singlecell_boxplot/*")


d5=data.frame()

for(i in seq_along(t)){
  print(t[i])
  
  if(str_detect(t[i],"cell")){
    
    d=read.csv(t[i],header=T)
  }
  else{
    print(t[i])
    d=read.csv(t[i],header=T)%>%
      filter(cancer=="yes")%>%
      select(program,data,Gene,TPM)
  }
  d5=rbind(d5,d)
}


#####################################other format
d5_aver=d5%>%
  filter(Gene!='SMAD9')%>%
  filter(Gene!='ID4')%>%
  group_by(program,data)%>%
  summarise(meanTPM=mean(TPM),n=n())

p=ggplot(d5_aver,aes(x=data,y=log2(meanTPM+1)))+
  geom_boxplot(fill="#e76f51",alpha=0.8)+
  geom_jitter()+
  theme_classic(base_size = 15)+
  ylab("Averaged exp for ID1,ID2 and ID3 Log2(TPM+1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

p1=ggplot(d5[d5$Gene=="ID1",],aes(x=data,y=log2(TPM+1)))+
  geom_boxplot(fill="#f4a261",alpha=0.5)+
  geom_jitter()+
  theme_classic(base_size = 15)+
  ylab("ID1 exp Log2(TPM+1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1

p2=ggplot(d5[d5$Gene=="ID2",],aes(x=data,y=log2(TPM+1)))+
  geom_boxplot(fill="#e9c46a",alpha=0.5)+
  geom_jitter()+
  theme_classic(base_size = 15)+
  ylab("ID2 exp Log2(TPM+1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2

p3=ggplot(d5[d5$Gene=="ID3",],aes(x=data,y=log2(TPM+1)))+
  geom_boxplot(fill="#2a9d8f",alpha=0.5)+
  geom_jitter()+
  theme_classic(base_size = 15)+
  ylab("ID3 exp Log2(TPM+1)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3

pdf("boxplot_datasets.pdf", width=12, height=12)

plot_grid(p,p1,p2,p3,nrow=2,rel_widths = c(2,2))

dev.off()