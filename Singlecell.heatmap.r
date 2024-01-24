#setwd("C:/Users/yzhang24/OneDrive - St. Jude Children's Research Hospital/Documents/sjproject/single_RNA")
setwd("Z:/ResearchHome/Groups/geelegrp/home/yzhang24/1_RA_BMP/3_RA_SMAD9_Project/5_single_RNA/")
library("ggplot2")
library(tidyverse)
library("DESeq2")
library(dplyr)
library("clusterProfiler")
library("enrichplot")
library(cowplot)
###############################################################################################Usage Description
##this script is used to extract interested gene expression from the programs annotations from single cell rnaseq data of 5 database
##and combined them together to explore how data looks like.
###############################################################################################

###############reading the TPM expression tables from bone marrow samples and primary tumor samples.
##1. reading tables from GOSH primary samples
###defined a list of genes that I am interested in

genelist=c("ID1","ID2","ID3","ID4","SMAD9")

t_gosh <- Sys.glob("./GOSH/GEP_*.csv")
t_dong <- Sys.glob("./dong/GEP_*.csv")
t_bmr <- Sys.glob("./BMR/reports_baseNMF_combined/GEP_*.csv")
t_human_pdx <- Sys.glob("./PDX/PDX_Human/GEP_*.csv")

d_gosh=data.frame()
d_dong=data.frame()
d_bmr=data.frame()
d_human_pdx =data.frame()

####################extract genes TPM from all programs from GOSH data set
for(i in seq_along(t_gosh)){
  
  d_temp=data.frame()
####integrate all tables into a data frame, is that possible?
  d=read.csv(t_gosh[i],header=T)
  row=0
  for (r in seq_along(d$Gene)){
    
    for (j in seq_along(genelist)) {
      
      if(d$Gene[r]==genelist[j]){
        #print(d[r,] )
        d_temp=rbind(d_temp,d[r,])
      }
    }
  }
  
  d_temp$table=basename(t_gosh[i])
  d_gosh=rbind(d_gosh,d_temp)

}
##################for dong dataset

for(i in seq_along(t_dong)){
  
  d_temp=data.frame()
  ####integrate all tables into a data frame, is that possible?
  d=read.csv(t_dong[i],header=T)
  row=0
  for (r in seq_along(d$Gene)){
    
    for (j in seq_along(genelist)) {
      
      if(d$Gene[r]==genelist[j]){
        #print(d[r,] )
        d_temp=rbind(d_temp,d[r,])
      }
    }
  }
  
  d_temp$table=basename(t_dong[i])
  d_dong=rbind(d_dong,d_temp)
}


####################extract genes TPM from all programs from BMR data set

for(i in seq_along(t_bmr)){
  
  d_temp=data.frame()
  ####integrate all tables into a data frame, is that possible?
  d=read.csv(t_bmr[i],header=T)
  row=0
  for (r in seq_along(d$Gene)){
    
    for (j in seq_along(genelist)) {
      
      if(d$Gene[r]==genelist[j]){
        #print(d[r,] )
        d_temp=rbind(d_temp,d[r,])
      }
    }
  }
  
  d_temp$table=basename(t_bmr[i])
  d_bmr=rbind(d_bmr,d_temp)
}

###########################extract genes TPM from all programs from human PDX data set

for(i in seq_along(t_human_pdx)){
  
  d_temp=data.frame()
  
  ####integrate all tables into a data frame, is that possible?
  d=read.csv(t_human_pdx[i],header=T)
  row=0
  for (r in seq_along(d$Gene)){
    
    for (j in seq_along(genelist)) {
      
      if(d$Gene[r]==genelist[j]){
        #print(d[r,] )
        d_temp=rbind(d_temp,d[r,])
      }
    }
  }
  
  d_temp$table=basename(t_human_pdx[i])
  d_human_pdx=rbind(d_human_pdx,d_temp)
}

#############################################Annotation and result plot output

a_gosh=read.csv("./GOSH_annotation.csv",header=T)
a_dong=read.csv("./dong_annotation.csv",header=T)
a_bmr=read.csv("./bone_marrow_annotation.csv",header = T)
a_human_pdx=read.csv("./humanPDX_annotation.csv",header = T)

a_gosh$program=as.character(a_gosh$program)
a_dong$program=as.character(a_dong$program)
a_bmr$program=as.character(a_bmr$program)
a_human_pdx$program=as.character(a_human_pdx$program)

d_gosh=d_gosh %>%
  mutate(program=as.character(str_match(d_gosh$table,"\\d+")))%>%
  mutate(data="GOSH")

d_dong=d_dong %>%
  mutate(program=as.character(str_match(d_dong$table,"\\d+")))%>%
  mutate(data="Dong")

d_bmr=d_bmr %>%
  mutate(program=as.character(str_match(d_bmr$table,"\\d+")))%>%
  mutate(data="BMR")

d_human_pdx=d_human_pdx %>%
  mutate(program=as.character(str_match(d_human_pdx$table,"\\d+")))%>%
  mutate(data="Human_PDX")


##################read processed cell line data & mouse tumor data TMP value in a format that can be combined with single cell data

d_gosh_anno=left_join(d_gosh,a_gosh,by="program")
d_dong_anno=left_join(d_dong,a_dong,by="program")
d_bmr_anno=left_join(d_bmr,a_bmr,by="program")
d_human_pdx_anno=left_join(d_human_pdx,a_human_pdx,by="program")

##################################################Output of target gene TPM
write.table(d_gosh_anno[,1:9],"./d_gosh_anno.IDs.TPM.table",quote = F,row.names = F,sep = "\t")
write.table(d_dong_anno[,1:9],"./d_dong_anno.IDs.TPM.table",quote = F,row.names = F,sep = "\t")
write.table(d_bmr_anno[,1:9],"./d_bmr_anno.IDs.TPM.table",quote = F,row.names = F,sep = "\t")
write.table(d_human_pdx_anno[,1:9],"./d_human_pdx_anno.IDs.TPM.table",quote = F,row.names = F,sep = "\t")

  #aa=read.table("./d_dong_anno.IDs.TPM.table",header = T,sep = "\t")


d_cell_anno=read.csv("./cell.mouse.format.csv",header=T)

d_bmr_anno=rbind(d_bmr_anno,d_cell_anno)%>%
            mutate(anno2=paste(program,annotation,sep="-"))%>%
            mutate(textcolor=ifelse(cancer=='yes',"red",ifelse(cancer=="no","black","blue")))

d_gosh_anno=rbind(d_gosh_anno,d_cell_anno) %>%
  mutate(anno2=paste(program,annotation,sep="-"))%>%
  mutate(textcolor=ifelse(cancer=='yes',"red",ifelse(cancer=="no","black","blue")))

d_dong_anno=rbind(d_dong_anno,d_cell_anno) %>%
  mutate(anno2=paste(program,annotation,sep="-"))%>%
  mutate(textcolor=ifelse(cancer=='yes',"red",ifelse(cancer=="no","black","blue")))

d_human_pdx_anno=rbind(d_human_pdx_anno,d_cell_anno) %>%
  mutate(anno2=paste(program,annotation,sep="-"))%>%
  mutate(textcolor=ifelse(cancer=='yes',"red",ifelse(cancer=="no","black","blue")))




#######################Heatmap for bone marrow
bmr_id1=d_bmr_anno[d_bmr_anno$Gene=='ID1',]
bmr_id2=d_bmr_anno[d_bmr_anno$Gene=='ID2',]
bmr_id3=d_bmr_anno[d_bmr_anno$Gene=='ID3',]

d_bmr_anno$anno2=factor(d_bmr_anno$anno2,levels =bmr_id1[order(bmr_id1$TPM,decreasing = F),'anno2'] )
id1color=bmr_id1[order(bmr_id1$TPM,decreasing = F),'textcolor']

p1=ggplot(d_bmr_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id1color))
p1

d_bmr_anno$anno2=factor(d_bmr_anno$anno2,levels =bmr_id2[order(bmr_id2$TPM,decreasing = F),'anno2'] )
id2color=bmr_id2[order(bmr_id2$TPM,decreasing = F),'textcolor']

p2=ggplot(d_bmr_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id2color))

p2

d_bmr_anno$anno2=factor(d_bmr_anno$anno2,levels =bmr_id3[order(bmr_id3$TPM,decreasing = F),'anno2'] )
id3color=bmr_id3[order(bmr_id3$TPM,decreasing = F),'textcolor']

p3=ggplot(d_bmr_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id3color))

p3

pdf("heatmap_bmr_allprogram.addcellline.id1,id2,id3.pdf", width=20, height=10)

plot_grid(p1,p2, p3,ncol=3)

dev.off()

########################Heatmap for GOSH

gosh_id1=d_gosh_anno[d_gosh_anno$Gene=='ID1',]
gosh_id2=d_gosh_anno[d_gosh_anno$Gene=='ID2',]
gosh_id3=d_gosh_anno[d_gosh_anno$Gene=='ID3',]

d_gosh_anno$anno2=factor(d_gosh_anno$anno2,levels =gosh_id1[order(gosh_id1$TPM,decreasing = F),'anno2'])
id1color=gosh_id1[order(gosh_id1$TPM,decreasing = F),'textcolor']

p1=ggplot(d_gosh_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id1color))
p1

d_gosh_anno$anno2=factor(d_gosh_anno$anno2,levels =gosh_id2[order(gosh_id2$TPM,decreasing = F),'anno2'])
id2color=gosh_id2[order(gosh_id2$TPM,decreasing = F),'textcolor']

p2=ggplot(d_gosh_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id2color))

p2


d_gosh_anno$anno2=factor(d_gosh_anno$anno2,levels =gosh_id3[order(gosh_id3$TPM,decreasing = F),'anno2'])
id3color=gosh_id3[order(gosh_id3$TPM,decreasing = F),'textcolor']

p3=ggplot(d_gosh_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id3color))

p3

pdf("heatmap_gosh_allprogram.addcellline.id1,id2,id3.pdf", width=20, height=10)

plot_grid(p1, p2,p3,ncol=3)

dev.off()
#################################################dong 
dong_id1=d_dong_anno[d_dong_anno$Gene=='ID1',]
dong_id2=d_dong_anno[d_dong_anno$Gene=='ID2',]
dong_id3=d_dong_anno[d_dong_anno$Gene=='ID3',]

d_dong_anno$anno2=factor(d_dong_anno$anno2,levels =dong_id1[order(dong_id1$TPM,decreasing = F),'anno2'])
id1color=dong_id1[order(dong_id1$TPM,decreasing = F),'textcolor']

p1=ggplot(d_dong_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id1color))

p1

d_dong_anno$anno2=factor(d_dong_anno$anno2,levels =dong_id2[order(dong_id2$TPM,decreasing = F),'anno2'])
id2color=dong_id2[order(dong_id2$TPM,decreasing = F),'textcolor']

p2=ggplot(d_dong_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id2color))

p2


d_dong_anno$anno2=factor(d_dong_anno$anno2,levels =dong_id3[order(dong_id3$TPM,decreasing = F),'anno2'])
id3color=dong_id3[order(dong_id3$TPM,decreasing = F),'textcolor']

p3=ggplot(d_dong_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id3color))

p3

pdf("heatmap_dong_allprogram.addcellline.id1,id2,id3.pdf", width=20, height=10)

plot_grid(p1, p2,p3,ncol=3)

dev.off()

###############################################PDX human
humanPDX_id1=d_human_pdx_anno[d_human_pdx_anno$Gene=='ID1',]


d_human_pdx_anno$anno2=factor(d_human_pdx_anno$anno2,levels =humanPDX_id1[order(humanPDX_id1$TPM,decreasing = F),'anno2'])
id1color=humanPDX_id1[order(humanPDX_id1$TPM,decreasing = F),'textcolor']

p1=ggplot(d_human_pdx_anno,aes(x=Gene,y=anno2,fill=log2(TPM+1)))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white",
                       high = "red", midpoint = 6,limits=c(0, 12))+
  theme_classic()+
  xlab("Log2(TPM+1)")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(colour = id1color))

p1

pdf("heatmap_humanPDX_allprogram.addcellline.id1,id2,id3.pdf", width=7, height=10)

plot_grid(p1,ncol=1)

dev.off()