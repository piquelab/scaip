#
library("rhdf5")
library(Seurat)
library(Matrix)
library(tidyverse)
library(cellranger)
library(parallel)
library(data.table)
#library(diem)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
theme_set(theme_grey())

rm(list=ls())
outdir <- "./3_compare_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


###################
### 1, expnames ###
###################
basefolder <- "/nfs/rprdata/julong/SCAIP/count/"
expNames <- dir(basefolder,"^SCAIP*")
folders <- paste0(basefolder, expNames, "/outs/filtered_feature_bc_matrix/", sep="")
ind <- dir.exists(folders) #ind <- file.info(folders)$isdir;ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
#cat(folders, sep="\n")
## remove bad emulsions
EXP_del <- c("SCAIP2-LPS-DEX", "SCAIP2-PHA-DEX", "SCAIP2-PHA-EtOH", "SCAIP3-LPS-DEX", "SCAIP3-PHA-DEX", "SCAIP3-PHA-EtOH")
folders <- folders[!expNames%in%EXP_del]
expNames <- names(folders)


a <- mclapply(expNames,function(ii){
  cat("#Loading ",ii,"\n")
  aux3 <- Read10X(folders[ii])
  colnames(aux3) <- gsub("-1","",colnames(aux3))
  return(aux3)
},mc.cores=5)
a <- do.call(cbind,a)
sc <- CreateSeuratObject(a,project="SCAIP")
##If not identical we may need to re-order columns. 
#stopifnot(identical(colnames(a),demuxlet$NEW_BARCODE))
###output
opfn <- "./3_compare_output/1_SCAIP.cellranger.rds"
write_rds(sc,opfn)



#####################################
### 1, compare reads not filtered ###
#####################################

######################################
### (1), compared Old vs New ###
######################################

rm(list=ls())
sc1 <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta1 <- sc1@meta.data 
New <- meta1%>%
       group_by(orig.ident)%>%
       summarise(ncell=n(),
                 reads=mean(nCount_RNA),
                 ngene=mean(nFeature_RNA))

sc2 <- read_rds("../SCAIP-ALL-2019.10.24/2and3_Kallisto_Diem_Output/1_SCAIP.kb.rds")
meta2 <- sc2@meta.data
meta2$NEW_BARCODE <- rownames(meta2)
Old <- meta2%>%
       group_by(orig.ident)%>%
       summarise(ncell=n(),
                 reads=mean(nCount_RNA),
                 ngene=mean(nFeature_RNA))
##ncell                       
comb <- Old%>%
        left_join(New, by="orig.ident")%>%
        mutate(BATCH=substring(orig.ident, 1, 6))
        

imin <- min(c(comb$ncell.x, comb$ncell.y))
imax <- max(c(comb$ncell.x, comb$ncell.y)) 
fig1 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Old")+ylab("New")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("#Number of Barcode")+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))
              
       
imin <- min(c(comb$reads.x, comb$reads.y))
imax <- max(c(comb$reads.x, comb$reads.y))       
fig2 <- ggplot(comb,aes(x=reads.x,y=reads.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Old")+ylab("New")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("#Reads per cell")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.8,0.25),
              legend.background=element_blank(), 
              plot.title=element_text(hjust=0.5))
              
png("./3_compare_output/FigureTmp1.1.png",width=1000, height=600, res=120)
plot_grid(fig1, fig2)
dev.off()

#######################################################
### (2), comparing kb+diem and cellranger(filtered) ###
#######################################################
rm(list=ls())
sc1 <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta1 <- sc1@meta.data 
kb <- meta1%>%group_by(orig.ident)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA))
             
sc2 <- read_rds("./2_comparison_output/1_SCAIP.cellranger.rds")
meta2 <- sc2@meta.data
Cellranger <- meta2%>%group_by(orig.ident)%>%
              summarise(ncell=n(),
                        reads=mean(nCount_RNA),
                        ngene=mean(nFeature_RNA))
                       
comb <- Cellranger%>%
        left_join(kb,by="orig.ident")%>%
        mutate(BATCH=substring(orig.ident, 1, 6))                                     
             
imin <- min(c(comb$ncell.x,comb$ncell.y))
imax <- max(c(comb$ncell.x,comb$ncell.y)) 
fig1 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Cellranger")+ylab("kb+DIEM")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("#Number of Barcode")+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))

### Reads per cell
imin <- min(c(comb$reads.x,comb$reads.y))
imax <- max(c(comb$reads.x,comb$reads.y)) 
fig2 <- ggplot(comb,aes(x=reads.x, y=reads.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Cellranger")+ylab("kb+DIEM")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        #guides(colour=guide_legend(ncol=3))+
        ggtitle("#Reads per cell")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.8,0.25),
              legend.background=element_blank(),
              plot.title=element_text(hjust=0.5))

png("./3_compare_output/FigureTmp1.2.png",width=1000, height=600, res=120)
plot_grid(fig1, fig2)
dev.off()


################################################################
### 2, after demuxlet to assign individual identity to cells ### 
###       compare matching and mismatching                   ###
################################################################


##############################################################################
### 2, Compared mismatching and matching barcode in Old data and New data  ###
##############################################################################
rm(list=ls())

### (1) compare msimatching and matching in old data
sc <- read_rds("../SCAIP-ALL-2019.10.24/2and3_Kallisto_Diem_Output/1_SCAIP.kb.rds")
meta <- sc@meta.data
meta$NEW_BARCODE <- rownames(meta)
#sample2 <- setNames(x1$sample2, x1$sample)
#batch2 <- setNames(x1$Group,x1$sample)

demuxlet <- read_rds("../SCAIP-ALL-2019.10.24/1_demuxlet_output/1_demuxlet_ALL.rds")
meta <- meta%>%inner_join(demuxlet,by="NEW_BARCODE")

x <- read.table("../CorrectBatch1-5.txt", header=T)
correctBatch0 <- paste(x$sample, x$Group, sep="_")
#sample2 <- setNames(x$sample2, x$sample)
batch2 <- setNames(x$Group,x$sample)

meta0 <- meta%>%
         mutate(comb=paste(meta$BEST.GUESS,meta$BATCH, sep="_"),
                matching=ifelse(comb%in%correctBatch0, 1, 0))
##0, mismatching and 1 is matching 
         
dd <- meta0%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH", "matching")
mydf0 <- dd%>%
         group_by(BEST.GUESS,matching)%>%
         summarise(ncell=n())%>%
         mutate(BATCH=batch2[BEST.GUESS])%>%
         group_by(BEST.GUESS)%>%
         mutate(Perc=ncell/sum(ncell))
         
fig0 <- ggplot(mydf0,aes(x=as.character(BEST.GUESS), y=ncell, fill=factor(matching)))+
        geom_bar(stat="identity",position="fill")+
        xlab("")+ylab("")+
        scale_fill_discrete("", labels=c("0"="mismatching", "1"="matching"))+
        facet_wrap(~factor(BATCH), ncol=2, scales="free")+
        #geom_text(aes(label=ncell),vjust=-0.5, size=3)+
        theme_bw()+
        theme(legend.position=c(0.85,0.2),
              legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=6))
                
png("./3_compare_output/FigureTmp2.1.png", width=700, height=700, res=120)
fig0
dev.off()


#########################################################
### (2), compare mismatching and matching in new data ###
#########################################################
rm(list=ls())
sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta <-sc@meta.data
demux <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
meta <- meta%>%inner_join(demux,by="NEW_BARCODE")

x <- read.table("../CorrectBatch1-6.txt", header=T)
correctBatch0 <- paste(x$sample, x$Group, sep="_")
#sample2 <- setNames(x$sample2, x$sample)
batch2 <- setNames(x$Group,x$sample)

meta0 <- meta%>%
         mutate(comb=paste(meta$BEST.GUESS,meta$BATCH, sep="_"),
                matching=ifelse(comb%in%correctBatch0, 1, 0))
                
dd <- meta0%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH", "matching")
mydf0 <- dd%>%
         group_by(BEST.GUESS,matching)%>%
         summarise(ncell=n())%>%
         mutate(BATCH=batch2[BEST.GUESS])%>%
         group_by(BEST.GUESS)%>%
         mutate(Perc=ncell/sum(ncell))
         
fig0 <- ggplot(mydf0,aes(x=as.character(BEST.GUESS), y=ncell, fill=factor(matching)))+
        geom_bar(stat="identity",position="fill")+
        xlab("")+ylab("")+
        scale_fill_discrete("", labels=c("0"="mismatching", "1"="matching"))+
        facet_wrap(~factor(BATCH), ncol=3, scales="free")+
        #geom_text(aes(label=ncell),vjust=-0.5, size=3)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=8))   
                
png("./3_compare_output/FigureTmp2.2.png", width=900, height=600, res=110)
fig0
dev.off()


#############################
### (3), compare matching ###
#############################
rm(list=ls())

###
### (A), Old data, all vs matching Barcodes
sc <- read_rds("../SCAIP-ALL-2019.10.24/2and3_Kallisto_Diem_Output/1_SCAIP.kb.rds")
meta <- sc@meta.data
meta$NEW_BARCODE <- rownames(meta)

demuxlet <- read_rds("../SCAIP-ALL-2019.10.24/1_demuxlet_output/1_demuxlet_ALL.rds")
meta <- meta%>%inner_join(demuxlet,by="NEW_BARCODE")

x <- read.table("../CorrectBatch1-5.txt", header=T)
correctBatch0 <- paste(x$sample, x$Group, sep="_")
#sample2 <- setNames(x$sample2, x$sample)
batch2 <- setNames(x$Group,x$sample)

meta0 <- meta%>%
         mutate(comb=paste(meta$BEST.GUESS,meta$BATCH, sep="_"),
                matching=ifelse(comb%in%correctBatch0, 1, 0))
                
allx <- meta0%>%
        group_by(orig.ident)%>%
        summarise(ncell=n(),
                  reads=mean(nCount_RNA))
          
matchx <- meta0%>%
          filter(matching==1)%>%
          group_by(orig.ident)%>%
          summarise(ncell=n(),
                    reads=mean(nCount_RNA))
          
comb <- matchx%>%
        inner_join(allx, by="orig.ident")%>%
        mutate(BATCH=substring(orig.ident,1,6))

imin <- min(c(comb$ncell.x, comb$ncell.y))
imax <- max(c(comb$ncell.x, comb$ncell.y)) 
fig1 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Matching")+ylab("All")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("Old data")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.75,0.3),
              legend.background=element_blank(),
              plot.title=element_text(hjust=0.5))        


### (B), New data, all vs matching Barcodes
sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta <-sc@meta.data
demux <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
meta <- meta%>%inner_join(demux,by="NEW_BARCODE")

x <- read.table("../CorrectBatch1-6.txt", header=T)
correctBatch0 <- paste(x$sample, x$Group, sep="_")
#sample2 <- setNames(x$sample2, x$sample)
batch2 <- setNames(x$Group,x$sample)

meta1 <- meta%>%
         mutate(comb=paste(meta$BEST.GUESS,meta$BATCH, sep="_"),
                matching=ifelse(comb%in%correctBatch0, 1, 0))

allx <- meta1%>%
        group_by(orig.ident)%>%
        summarise(ncell=n(),
                  reads=mean(nCount_RNA))
          
matchx <- meta1%>%
          filter(matching==1)%>%
          group_by(orig.ident)%>%
          summarise(ncell=n(),
                    reads=mean(nCount_RNA))
          
comb <- matchx%>%
        inner_join(allx, by="orig.ident")%>%
        mutate(BATCH=substring(orig.ident,1,6))

imin <- min(c(comb$ncell.x, comb$ncell.y))
imax <- max(c(comb$ncell.x, comb$ncell.y)) 
fig2 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Matching")+ylab("All")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("New data")+
        theme_bw()+
        theme(legend.title=element_blank(),
              legend.position=c(0.75,0.3),
              legend.background=element_blank(),
              plot.title=element_text(hjust=0.5))                   
                

### (C), old vs new for matching barcodes
Old <- meta0%>%
       filter(matching==1)%>%
       group_by(orig.ident)%>%
       summarise(ncell=n(),
                 reads=mean(nCount_RNA))
New <- meta1%>%
       filter(matching==1)%>%
       group_by(orig.ident)%>%
       summarise(ncell=n(),
                 reads=mean(nCount_RNA))
comb <- New%>%
        inner_join(Old, by="orig.ident")%>%
        mutate(BATCH=substring(orig.ident,1,6))
imin <- min(c(comb$ncell.x, comb$ncell.y))
imax <- max(c(comb$ncell.x, comb$ncell.y)) 
fig3 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("New")+ylab("Old")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("Old vs New(matching)")+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5))

png("./3_compare_output/FigureTmp2.3.png", width=1200, height=500, res=120)
plot_grid(fig1, fig2, fig3,ncol=3)
dev.off()
               

#########################################################
### 3, after demuxlet, summary barcodes by individual ### 
#########################################################
rm(list=ls())   


        
###
### 3, demuxlet, Number of SNPs
rm(list=ls())
sc1 <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta1 <- sc1@meta.data 
New <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
New <- meta1%>%inner_join(New,by="NEW_BARCODE")

fig1 <- ggplot(New,aes(x=NUM.SNPS))+
        geom_histogram(color="#e9ecef", alpha=0.6, position="identity")+
        xlab("NUM.SNPs")+theme_bw() 
fig2 <- ggplot(New,aes(x=NUM.READS))+
        geom_histogram(color="#e9ecef", alpha=0.6, position="identity")+
        xlab("NUM.Reads")+theme_bw() 
        
png("./3_compare_output/FigureTmp3.0.png", width=1000, height=500, res=120)
plot_grid(fig1, fig2,ncol=3)
dev.off()  

x <- read.table("../CorrectBatch1-6.txt", header=T)
correctBatch0 <- paste(x$sample, x$Group, sep="_")
#sample2 <- setNames(x$sample2, x$sample)
batch2 <- setNames(x$Group,x$sample)

### Old data
sc <- read_rds("../SCAIP-ALL-2019.10.24/2and3_Kallisto_Diem_Output/1_SCAIP.kb.rds")
meta <- sc@meta.data
meta$NEW_BARCODE <- rownames(meta)
demux <- read_rds("../SCAIP-ALL-2019.10.24/1_demuxlet_output/1_demuxlet_ALL.rds")
meta0 <- meta%>%inner_join(demux,by="NEW_BARCODE")
old0 <- meta0%>%
        group_by(BEST.GUESS)%>%
        summarise(ncell=n())
        
### New data 
sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta <-sc@meta.data
demux <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
meta1 <- meta%>%inner_join(demux,by="NEW_BARCODE")
new0 <- meta1%>%
        group_by(BEST.GUESS)%>%
        summarise(ncell=n())

comb <- old0%>%
        left_join(new0,by="BEST.GUESS")%>%
        mutate(BATCH=batch2[BEST.GUESS])
        
imin <- min(c(comb$ncell.x,comb$ncell.y))
imax <- max(c(comb$ncell.x,comb$ncell.y)) 
fig1 <- ggplot(comb,aes(x=ncell.x,y=ncell.y, colour=factor(BATCH)))+
        geom_point()+
        xlab("Old")+ylab("New")+xlim(imin,imax)+ylim(imin,imax)+
        geom_abline(intercept=0,slope=1,colour="blue")+
        ggtitle("#Number of Barcode")+
        theme_bw()+
        theme(legend.title=element_blank(),
              plot.title=element_text(hjust=0.5)) 
          
png("./3_compare_output/FigureTmp3.1.png", width=600, height=500, res=120)
fig1
dev.off()     
        
        
        
          


#############
### 3, compare demuxlet, New vs Old
##########
#rm(list=ls())
#Old <- read_rds("../SCAIP-ALL-2019.10.24/1_demuxlet_output/1_demuxlet_ALL.rds")
#New <- read_rds("./1_demux2_output/1_demux_filter0.ALL.rds")

#old0 <- Old%>%select("NEW_BARCODE","BEST.GUESS", "EXP", "BATCH")
#new0 <- New%>%select("NEW_BARCODE","BEST.GUESS", "EXP", "BATCH")
#comb <- old0%>%inner_join(new0,by="NEW_BARCODE")

### By experiments
#Old0 <- Old%>%
#        group_by(EXP)%>%
#        summarise(ncell=n(),
#                  reads=mean(NUM.READS))
#New0 <- New%>%
#        group_by(EXP)%>%
#        summarise(ncell=n(),
#                  reads=mean(NUM.READS))
#Comb <- Old0%>%
#        left_join(New0,by="EXP")
#
## number of cells        
#imin <- min(c(Comb$ncell.x,Comb$ncell.y))
#imax <- max(c(Comb$ncell.x,Comb$ncell.y)) 

#fig0 <- ggplot(Comb,aes(x=ncell.x,y=ncell.y, colour=factor(EXP)))+
#        geom_point()+
#        xlab("Old")+ylab("New")+
#        xlim(imin,imax)+ylim(imin,imax)+
#        guides(colour=guide_legend(ncol=3))+
#        geom_abline(intercept=0,slope=1,colour="blue")+
#        ggtitle("#Number of Barcode")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.text=element_text(size=6), 
#              plot.title=element_text(hjust=0.5))
#png("./2_comparison_output/FigureTmp3.1.png",width=1000, height=600, res=120)
#fig0
#dev.off() 

## reads per cell barcode
#imin <- min(c(Comb$reads.x,Comb$reads.y))
#imax <- max(c(Comb$reads.x,Comb$reads.y))
#fig0 <- ggplot(Comb,aes(x=reads.x,y=reads.y, colour=factor(EXP)))+
#        geom_point()+
#        xlab("Old")+ylab("New")+
#        xlim(imin, imax)+ylim(imin,imax)+
#        guides(colour=guide_legend(ncol=3))+
#        geom_abline(intercept=0,slope=1,colour="blue")+
#        ggtitle("#Reads per cell")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.text=element_text(size=6), 
#              plot.title=element_text(hjust=0.5))
#png("./2_comparison_output/FigureTmp3.2.png",width=1000, height=600, res=120)
#fig0
#dev.off()

### by individuals
#Old0 <- Old%>%
#        group_by(BEST.GUESS)%>%
#        summarise(ncell=n(),
#                  reads=mean(NUM.READS))
#New0 <- New%>%
#        group_by(BEST.GUESS)%>%
#        summarise(ncell=n(),
#                  reads=mean(NUM.READS))
#Comb <- Old0%>%
#        left_join(New0,by="BEST.GUESS")
#
###        
#imin <- log10(min(c(Comb$ncell.x, Comb$ncell.y)))
#imax <- log10(max(c(Comb$ncell.x, Comb$ncell.y)))                
#fig0 <- ggplot(Comb,aes(x=log10(ncell.x), y=log10(ncell.y), colour=factor(BEST.GUESS)))+
#        geom_point()+
#        xlab("Old")+ylab("New")+
#        xlim(imin,imax)+ylim(imin,imax)+
#        guides(colour=guide_legend(ncol=5))+
#        geom_abline(intercept=0,slope=1,colour="blue")+
#        ggtitle("#Number of Barcode")+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              legend.text=element_text(size=6), 
#              plot.title=element_text(hjust=0.5))
#png("./2_comparison_output/FigureTmp3.3.png",width=1000, height=600, res=120)
#fig0
#dev.off()            

