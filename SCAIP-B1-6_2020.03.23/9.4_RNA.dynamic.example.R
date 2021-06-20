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
    yi <- sum(y[s0:s1])
    zi <- mean(z[s0:s1])  
    d <- data.frame(z=zi, y=yi)
    cvt2 <- rbind(cvt2, d)  
    s0 <- s0+step
  }
  cvt2
}



###
### required constant
option <- "DiagLDA2"

grch38_unq <- grch38%>%filter(grepl("protein_coding", biotype))%>%
    distinct(ensgene,.keep_all=TRUE)%>%
    dplyr::select(ensgene, symbol)
#res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")

contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"), "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))

LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")


###
### Read counts data
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
X <- sc@assays$RNA@data
X <- X[grepl("S-",rownames(X)),]
rn <- gsub("S-|\\..*", "", rownames(X))
rownames(X) <- rn
meta <- sc@meta.data%>%dplyr::select(NEW_BARCODE,BEST.GUESS)


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


### read LDA-interaction eQTL
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/DiagLDA2/lm_results2/signif_pairs-SNP-eGene.txt"
signif <- read.table(fn,header=T)
## pick up SNP-gene pairs
signif2 <-signif%>%filter(symbol=="IFNG")
for (k in 1:nrow(signif2)) {
k <- 1
oneMCl <- signif2$MCls[k]
lda <- signif2$LDA[k]
i_lda <- which(LDA==lda)
contrast <- contrast_ls[[lda]]
gene0 <- signif2$ENSG2[k]
symbol0 <- signif2$symbol[k]
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
   mutate(z=meta2[,i_lda],y=X[gene0,NEW_BARCODE])

sample <- unique(cvt$BEST.GUESS)
cvt2 <- map_dfr(contrast, function(oneTreat){
  tmp <- map_dfr(sample, function(ii){
    #cat(oneTreat, ii, "\n")  
    cvt0 <- cvt%>%filter(treats==oneTreat, BEST.GUESS==ii)
    cvt0 <- slideFun(cvt0, win=0.2, step=0.05)
    cvt0 <- as.data.frame(cvt0)%>%mutate(treats=oneTreat, BEST.GUESS=ii, dosage=gen[ii])
    cvt0
  })
  tmp
})

cvt2$y <- cvt2$y/max(cvt2$y)

###
### (2) plot curve
ref <- strsplit(varID, ":")[[1]][3]
alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])

cvt0 <- cvt2%>%drop_na(y)%>%filter(treats==contrast[2])
fig2 <- ggplot(cvt0, aes(x=z, y=y, colour=factor(dosage)))+
  geom_point(size=0.3)+  
  geom_smooth(method="loess", size=0.3, se=F)+
  scale_colour_manual("",
    values= c("0"="#4daf4a", "1"="#984ea3","2"="#ff7f00"),
    labels=c("0"=paste(ref, "/", ref, sep=""),
             "1"=paste(ref, "/", alt, sep=""),
             "2"=paste(alt, "/", alt, sep="")))+  
  scale_y_continuous("Relative changes", limits=c(0,max(cvt0$y)))+
  xlab(paste("LDA_", lda, sep=""))+
  ggtitle(bquote(~italic(.(symbol0))~" in "~.(oneMCl)~"_"~.(contrast[2]) ))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5, size=10),
        axis.title=element_text(size=8))
##
outdir <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
                "/Example_gene.gr.dosage/", sep="")
if (!file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
    
figfn <- paste(outdir, "Figure1_", oneMCl, ".lda", lda, ".trt",
   contrast[2], "_", symbol0, ".", varID, ".fitting.png", sep="")
png(figfn, width=400, height=400, res=120)
print(fig2)
dev.off()

}###
## covid_list <- c("OAS3", "TAC4", "DPP9", "RAVER1", "IFNAR2", "THBS3",
##                 "SCN1A", "LZTFL1", "FOXP4", "TMEM65", "ABO", "OAS1",
##                 "KANSL1", "RPL24", "DNAH5", "PLEKHA4")
