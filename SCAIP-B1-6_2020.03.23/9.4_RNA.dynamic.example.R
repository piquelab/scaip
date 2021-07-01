##
###
library("rhdf5")
library("corpcor")
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
library(purrr)
library(furrr)
##
library(Seurat)
#library(SeuratDisk)
#library(harmony)
library(annotables) 
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(viridis)
theme_set(theme_grey())

rm(list=ls())


###
###
slideFun <- function(cvt, win=0.1, step=0.001){
###
  cvt <- cvt%>%arrange(z)
    
  win <- trunc(nrow(cvt)*win)
  step <- trunc(nrow(cvt)*step)

  nlen <- nrow(cvt) 
  z <- cvt$z
  y <- cvt$y
  cvt2 <- NULL
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    d <- cvt[s0:s1,]%>%group_by(dosage)%>%
       summarise(y=mean(y),.groups="drop")%>%as.data.frame()   
    zi <- mean(z[s0:s1])  
    d$z_ave <- zi
    cvt2 <- rbind(cvt2, d)  
    s0 <- s0+step
  }
  #cvt <- cbind(cvt, cvt2)
  cvt2
}


###
slideFun2 <- function(cvt, win=0.1, step=0.001){
###
  cvt <- cvt%>%arrange(z)
    
  win <- trunc(nrow(cvt)*win)
  step <- trunc(nrow(cvt)*step)

  nlen <- nrow(cvt) 
  z <- cvt$z
  y <- cvt$y
  cvt2 <- NULL
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    d <- cvt[s0:s1,]%>%group_by(dosage)%>%
       summarise(y=mean(y),.groups="drop")%>%as.data.frame()   
    zi <- mean(z[s0:s1])  
    d$z_ave <- zi
    cvt2 <- rbind(cvt2, d)  
    s0 <- s0+step
  }
  #cvt <- cbind(cvt, cvt2)
  cvt2
}


###
### required constant
## grch38_unq <- grch38%>%filter(grepl("protein_coding", biotype))%>%
##     distinct(ensgene,.keep_all=TRUE)%>%
##     dplyr::select(ensgene, symbol)
## #res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")

option <- "DiagLDA2"
###
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"), "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
###
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
treats <- c("LPS-EtOH", "LPS-DEX", "PHA-EtOH", "PHA-DEX")
dataset <- data.frame(MCls=rep(MCls, each=4),
   LDA=rep(LDA, times=4), treats=rep(treats, times=4))
col1 <- c("CTRL"="#828282",
   "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")

sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
X <- sc@assays$RNA@data
X <- X[grepl("S-",rownames(X)),]
rn <- gsub("S-|\\..*", "", rownames(X))
rownames(X) <- rn
meta <- sc@meta.data%>%dplyr::select(NEW_BARCODE, BEST.GUESS)


###
### read dosage 
genfn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/SCAIP1-6_filtered.vcf.txt.gz"
dosages <- fread(genfn, data.table=F, stringsAsFactors=F)
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/sample.txt"
sample <- read.table(fn)
names(dosages) <- c("chr", "pos", "varID", sample$V1)
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID
#dosages <- dosages[!duplicated(rownames(dosages)),]


## oneMCl <- "Tcell"
## lda <- "LPS"
## i_lda <- which(LDA==lda)
## contrast <- contrast_ls[[lda]]

##

i <- 16
oneMCl <- dataset[i,1]
lda <- dataset[i,2]
treat <- dataset[i,3]
contrast <- contrast_ls[[lda]]
treat0 <- contrast[1]
treat1 <- contrast[2]
ii <- which(LDA==lda)
col0 <- col1[treat]

### read LDA-interaction eQTL
fn <- paste("/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/DiagLDA2/",
   "plots/eGene_SNP.", i, "_", oneMCl, ".lda", lda, ".trt", treat, ".txt", sep="")
signif <- read.table(fn, header=T)
## pick up SNP-gene pairs
signif2 <-signif%>%filter(symbol=="KLRC1")
for (k in 1:nrow(signif2)) {

   ENSG <- signif2$ENSG[k] 
   ENSG2 <- signif2$ENSG2[k]
   symbol <- signif2$symbol[k]
   varID <- signif2$varID[k]
   gen <- dosages[varID,-(1:3)]
   gen <- round(as.numeric(gen[1,]))
   names(gen) <- names(dosages)[-(1:3)]

###
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
      "/1_MCls.", oneMCl, ".DLDA.rds", sep="")
   meta2 <- read_rds(fn)%>%left_join(meta, by="NEW_BARCODE")

## meta2$dosage <- gen[meta2$BEST.GUESS]
## x <- meta2%>%filter(treats==contrast[2])%>%
##   distinct(BEST.GUESS,.keep_all=TRUE) #%>%

  cvt <- meta2%>%dplyr::select(NEW_BARCODE, treats, BEST.GUESS)%>%
     mutate(z=meta2[,ii], y=X[ENSG2,NEW_BARCODE],
            dosage=gen[meta2$BEST.GUESS])
###
  ## genotype <- sort(unique(cvt$dosage))  
  cvt2 <- map_dfr(contrast, function(oneTreat){
    ## tmp <- map_dfr(genotype, function(ii){
    ## cat(oneTreat, ii, "\n")  
     cvt0 <- cvt%>%filter(treats==oneTreat)
     cvt0 <- slideFun(cvt0, win=0.1, step=0.01)
     cvt0 <- as.data.frame(cvt0)%>%mutate(treats=oneTreat)
     cvt0
    ## })
    ## tmp
  })
    
  #cvt2$y <- cvt2$y/max(cvt2$y)

### (2) plot curve
   ref <- strsplit(varID, ":")[[1]][3]
   alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])

   cvt0 <- cvt2%>%drop_na(y)%>%filter(treats==contrast[2])
   cvt0$y2 <- cvt0$y/max(cvt0$y) 
   fig2 <- ggplot(cvt0, aes(x=z_ave, y=y2))+
     geom_point(aes(shape=factor(dosage)), colour=col0, size=1)+  
     geom_smooth(aes(linetype=factor(dosage)), method="loess",
        colour=col0, size=0.5, se=F)+
     scale_linetype_manual("genotype",
        values= c("0"="dotted", "1"="dotdash", "2"="solid"),
        labels=c("0"=paste(ref, "/", ref, sep=""),
                 "1"=paste(ref, "/", alt, sep=""),
                 "2"=paste(alt, "/", alt, sep="")))+
     scale_shape_manual("genotype",
        values= c("0"=9, "1"=8, "2"=25),
        labels=c("0"=paste(ref, "/", ref, sep=""),
                 "1"=paste(ref, "/", alt, sep=""),
                 "2"=paste(alt, "/", alt, sep="")))+
     scale_y_continuous("Relative changes", limits=c(0,max(cvt0$y2)))+
     xlab(paste("LDA_", lda, sep=""))+
     ggtitle(bquote(~italic(.(symbol))~" in "~.(oneMCl)~"_"~.(contrast[2]) ))+
     theme_bw()+
     theme(plot.title=element_text(hjust=0.5, size=10),
         axis.title=element_text(size=8))
##
outdir <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
                "/Example_gene.gr.dosage/", sep="")
if (!file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
    
figfn <- paste(outdir, oneMCl, ".lda", lda, ".trt",
   contrast[2], "_", ENSG, "_", symbol, "_", varID, ".2.fitting.png", sep="")
png(figfn, width=400, height=400, res=120)
print(fig2)
dev.off()

}

###
## covid_list <- c("OAS3", "TAC4", "DPP9", "RAVER1", "IFNAR2", "THBS3",
##                 "SCN1A", "LZTFL1", "FOXP4", "TMEM65", "ABO", "OAS1",
##                 "KANSL1", "RPL24", "DNAH5", "PLEKHA4")
