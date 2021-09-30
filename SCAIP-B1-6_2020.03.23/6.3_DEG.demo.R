###
library(Matrix)
## library(MASS)
## library(scales)
library(data.table)
library(tidyverse)
library(Seurat)
##
library(annotables)
library(biobroom)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(viridis)
##

rm(list=ls())
theme_set(theme_grey())
outdir <- "./6_DEG.CelltypeNew_output/Filter2_pub/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)



#######################
#### demo gene show ###
#######################

###
### count data
## sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
## count <- sc@assays$RNA@counts
## meta <- sc@meta.data

## grch38_unq <- grch38%>%
v##    distinct(ensgene,.keep_all=T)%>%
##    dplyr::select(ensgene, symbol, chr, biotype)
## vars <- data.frame(rn=rownames(count))%>%
##    mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn),
##           ensgene2=gsub("[SU]-", "", rn),
##           uns=grepl("S-",rn), rnz=rowSums(count))%>%
##           left_join(grch38_unq, by="ensgene")

## autosome <- as.character(1:22)
## varsSel <- vars%>%filter(uns, rnz>20, chr%in%autosome)#, grepl("protein_coding", biotype))
## ### Out of 17,697 protein coding gene, 15,770 with reads>20
## ### Out of 42,554 genes, 30,972 whole genes>20
## X <- count[varsSel$rn,]
## rownames(X) <- varsSel$ensgene
## size <- colSums(X)/1e+06

## ### obs, cells information
## obs <- meta%>%filter(MCls=="Tcell", grepl("PHA", treats))%>%
##    mutate(bti=paste(MCls, treats, BEST.GUESS, BATCH, sep="_"))%>%
##    dplyr::select(NEW_BARCODE, bti, MCls, treats)
## dd <- obs%>%group_by(bti)%>%summarise(ncell=n(), .groups="drop")%>%
##    separate(col=bti, c("MCls", "treats", "BEST.GUESS", "BATCH"),
##             remove=FALSE, sep="_")%>%
##    arrange(desc(ncell))


## ######################
## ### candidate gene ###
## ######################
## ###
## resDEG <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
## resDEG2 <- resDEG%>%drop_na(qval)%>%
##    filter(MCls=="Tcell", contrast=="PHA-DEX")%>%
##    dplyr::select(rn, baseMeanHat) 
## x <- resDEG%>%drop_na(qval)%>%
##    filter(MCls=="Tcell", contrast=="PHA-DEX", abs(beta)>0.5, qval<0.1)%>%
##    dplyr::select(rn, beta, qval, baseMeanHat)%>%
##    arrange(desc(abs(beta)))

## ##
## resDMG <- read.table("./10_RNA.Variance_output/tmp9/4_mu.meta", header=TRUE)
## resDMG2 <- resDMG%>%drop_na(qval)%>%
##    filter(MCls=="Tcell", contrast=="PHA-DEX")%>%
##    dplyr::select(rn, gene, beta, pval, qval)

## resDVG <- read.table("./10_RNA.Variance_output/tmp9/3_phiNew.meta", header=TRUE)
## resDVG2 <- resDVG%>%drop_na(qval)%>%
##    filter(MCls=="Tcell", contrast=="PHA-DEX")%>%
##    dplyr::select(rn,  beta, pval, qval)

## res2 <- resDMG2%>%inner_join(resDVG2, by="rn")
## qx <- res2$qval.x
## qy <- res2$qval.y
## bx <- abs(res2$beta.x)
## by <- abs(res2$beta.y)
## gr <- rep(0, nrow(res2))
## gr[qx<0.1&bx>0.5] <- 1
## gr[qy<0.1&by>0.5] <- 2
## gr[(qx<0.1&bx>0.5)&(qy<0.1&by>0.5)] <- 3
## res2$gr <- gr

## res3 <- res2%>%filter(gr==1)%>%left_join(resDEG2, by="rn")%>%
##    arrange(desc(baseMeanHat))

## res3[1:20,c("gene","beta.x", "beta.y", "baseMeanHat")]

## ###
## ###
## gene0 <- "ENSG00000166710"
## ## btiSel <- c("Tcell_PHA-EtOH_AL-249_SCAIP5", "Tcell_PHA-DEX_AL-249_SCAIP5")
## ## df2 <- map_dfr(btiSel, function(ii){
## ##    treats <- gsub("Tcell_|_AL.*", "", ii) 
## ##    celli <- obs%>%filter(bti==ii)%>%dplyr::pull(NEW_BARCODE)
## ##    xi <- X[gene0, celli]
## ##    xi <- xi#/size[celli]
## ##    dfi <- data.frame(y=xi, treats=treats)
## ##    dfi
## ## })
## btiSel <- dd%>%filter(BATCH=="SCAIP6")%>%dplyr::pull(bti)
## obs2 <- obs%>%filter(bti%in%btiSel)%>%dplyr::select(NEW_BARCODE,treats)
## xi <- X[gene0, obs2$NEW_BARCODE]
## xi <- xi#/size[celli]
## dfi <- data.frame(y=xi, treats=obs2$treats)


## col1 <- c("CTRL"="#828282",
##   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

## df2 <- dfi%>%filter((y>6)&(y<100))
## p2 <- ggplot(df2, aes(x=factor(treats),y=y, fill=treats))+
##    geom_violin()+
##    scale_fill_manual(values=c("PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4"))+ 
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank())

## figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure7.png"
## png(figfn, width=480, height=480, res=120)
## print(p2)
## dev.off()


#######################
### simulation data ###
#######################

### differential expression
x <- seq(-4, 4, length=100)
y <- dnorm(x)
d1 <- data.frame(x=x, y=y)
###
x <- x-2
y <- dnorm(x, -2, 1)
d2 <- data.frame(x=x, y=y)

df2 <- rbind(d1, d2)
df2$contrast <- rep(c("PHA", "PHA-DEX"), each=100)

p1 <- ggplot(df2)+
   geom_line(aes(x=x, y=y, colour=factor(contrast)), size=1)+
   scale_color_manual(values=c("PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
   ggtitle("Shift in expression")+
   ylab("Density")+ylim(0,0.45)+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text=element_blank(),
         axis.ticks=element_blank(),
         plot.title=element_text(hjust=0.5, size=12))

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.3_DEG.png"
png(figfn, width=325, height=250, res=120)
print(p1)
dev.off()


### differential dispersion
x <- seq(-4, 4, leng=100)
y <- dnorm(x, 0, 1)
d1 <- data.frame(x=x, y=y)
###
x <- x
y <- dnorm(x, 0, 2)
d2 <- data.frame(x=x, y=y)

df2 <- rbind(d1, d2)
df2$contrast <- rep(c("PHA", "PHA-DEX"), each=100)

p2 <- ggplot(df2)+
   geom_line(aes(x=x, y=y, colour=factor(contrast)), size=1)+
   scale_color_manual(values=c("PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
   ggtitle("Shift in variability")+
   ylab("Density")+ 
   theme_bw()+
   theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text=element_blank(),
         axis.ticks=element_blank(),
         plot.title=element_text(hjust=0.5, size=12))

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.4_DVG.png"
png(figfn, width=325, height=250, res=120)
print(p2)
dev.off()


###
###
x <- seq(-4, 4, leng=100)
y <- dnorm(x, 0, 1)
d1 <- data.frame(x=x, y=y)
###
x <- seq(-6,6,length=100)
y <- dnorm(x, 0, 2)
d2 <- data.frame(x=x, y=y)

df2 <- rbind(d1, d2)
df2$contrast <- rep(c("PHA", "PHA-DEX"), each=100)

p2 <- ggplot(df2)+
   geom_line(aes(x=x, y=y, colour=factor(contrast)), size=1)+
   scale_color_manual(values=c("PHA"="#FFB9B9", "PHA-DEX"="#831717"))+
   ## ggtitle("Equal mean expression")+
   theme_void()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5,size=12)) 

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.4.2_DVG.png"
png(figfn, width=350, height=250, res=120)
print(p2)
dev.off()


###
###
x <- seq(-8, 8, length=300)
y <- dnorm(x, 0, 2)
d1 <- data.frame(x=x, y=y)

p1 <- ggplot(d1)+
   geom_line(aes(x=x, y=y), colour="#828282", size=1)+
   theme_void()

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.4.3._CTRL.png"
png(figfn, width=400, height=200, res=120)
print(p1)
dev.off()

###
x <- seq(-4,4,length=100)
y <- dnorm(x, 0, 0.8)
d2 <- data.frame(x=x, y=y)
p2 <- ggplot(d2)+
   geom_line(aes(x=x, y=y), colour="#fb9a99", size=1)+
   theme_void()

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.4.4._LPS.png"
png(figfn, width=280, height=250, res=120)
print(p2)
dev.off()

#################
### LDA plots ###
################

fn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA2/1_MCls.Tcell.DLDA.rds"
df <- read_rds(fn)

df2 <- df%>%filter(treat2%in%c("PHA", "PHA-DEX"))%>%
   dplyr::select("z_PHA", "z_PHA-DEX", "treat2")%>%
   dplyr::rename("x"="z_PHA-DEX", "y"="z_PHA")

p1 <- ggplot(df2)+
   geom_point(aes(x=x, y=y, colour=factor(treat2)), size=0.1)+
   ## geom_density(data=df2,aes(x, colour=factor(treat2)), size=1, fill=NA)+
   scale_colour_manual(values=c("CTRL"="#828282",
                                "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
   labs(x="LDA of PHA+DEX", y="LDA of PHA")+
   theme_bw()+
   theme(legend.position="none")
p2 <- ggMarginal(p1, groupColour=TRUE, groupFill=TRUE, size=3)

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.5_DLDA.png"
png(figfn, width=480, height=480, res=120)
print(p2)
dev.off()



############################
### An example of reQTL  ###
############################

rm(list=ls())
### treatment 1
prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/"
fn <- paste(prefix, "eQTL/normalized_GE_residuals/for-plotting/NKcell_PHA-EtOH_3GEPCs-rmed.bed", sep="")
X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
GE_treat1 <- X[,5:ncol(X)]
rownames(GE_treat1) <- X$ID


### treatment 2
fn <- paste(prefix, "eQTL/normalized_GE_residuals/for-plotting/NKcell_PHA-DEX_3GEPCs-rmed.bed", sep="")
X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
GE_treat2 <- X[,5:ncol(X)]
rownames(GE_treat2) <-  X$ID


### dosage
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/SCAIP1-6_filtered.vcf.txt.gz"
dosages <- fread(fn, header=F, data.table=F, stringsAsFactors=FALSE)
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/sample.txt"
IDs <- read.table(fn,stringsAsFactors=FALSE)[,1]
colnames(dosages) <- c("chr", "pos", "varID", IDs)
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID

###
rn <- rownames(GE_treat1)
x <- grch38%>%filter(symbol=="CAPZB")
ENSG <- rn[grepl(x[1,1], rn)]
varID <- dosages[grepl("rs6696220", dosages$varID), "varID"]    
symbol <- "CAPZB"


#################
### plot data ###
#################

### treat 1, PHA-EtOH
y <- GE_treat1[ENSG,]
x <- dosages[varID, names(y)]
d1 <- data.frame(y=as.numeric(y), x=round(as.numeric(x)), treats="PHA")
### treat 2, PHA-DEX
y <- GE_treat2[ENSG,]
x <- dosages[varID, names(y)]
d2 <- data.frame(y=as.numeric(y), x=round(as.numeric(x)), treats="PHA-DEX")

df <- rbind(d1, d2)

###
ref <- strsplit(varID, ":")[[1]][3]
alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])
p1 <- ggplot(data=df, aes(x=factor(x), y=y))+
   geom_boxplot(aes(fill=treats), outlier.shape=NA)+
   scale_x_discrete("rs6696220", labels=c("0"=paste(ref, "/", ref, sep=""),
      "1"=paste(ref, "/", alt, sep=""),
      "2"=paste(alt, "/", alt, sep="")))+
   scale_fill_manual("",
      values=c("PHA"="#a6cee3", "PHA-DEX"="#1f78b4"), guide="none")+
   ## geom_smooth(data=df, aes(x=as.numeric(x), y=as.numeric(y)),
   ##    method='lm', formula=y~as.numeric(x), se=FALSE)+
   geom_jitter(width=0.25, size=1)+
   facet_grid(.~treats,
      labeller=labeller(treats=c("PHA"="PHA", "PHA-DEX"="PHA+DEX")) )+
   ylab(bquote(~italic(.(symbol))~" normalized gene expresssion"))+
   ## ggtitle("NK cell")+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=10),
         strip.text.x=element_text(size=12))## ,
         ## plot.title=element_text(hjust=0.5, size=12))
figfn <- paste("./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.6_",
   symbol, "-", varID, ".png", sep="")
png(figfn, width=500, height=500, res=120)
print(p1)
dev.off()

write_rds(p1, "./6_DEG.CelltypeNew_output/Filter2_pub/1.6_reQTL.rds")



##########################################
### An example of vQTLs while not eQTL ### 
##########################################

prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/"

### dosage
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/SCAIP1-6_filtered.vcf.txt.gz"
dosages <- fread(fn, header=F, data.table=F, stringsAsFactors=FALSE)
fn <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/check_output/sample.txt"
IDs <- read.table(fn, stringsAsFactors=FALSE)[,1]
colnames(dosages) <- c("chr", "pos", "varID", IDs)
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID

##################
### expression ###
##################

fn <- paste(prefix, "eQTL-on-mean/qnormed_mean/Tcell_PHA-DEX_mean.bed.gz", sep="")
X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
GE <- X[,5:ncol(X)]
rownames(GE) <- X$ID

##convariates
fn <- paste(prefix, "eQTL-on-mean/covariates/Tcell_PHA-DEX_covs_8PCs.txt", sep="")
cv.mean <- read.table(fn, header=TRUE, row.names=1)
names(cv.mean) <- gsub("\\.", "-", names(cv.mean))


#######################
### dispersion data ###
#######################

fn <- paste(prefix, "dispersionQTL/qnormed_dispersion/Tcell_PHA-DEX_dispersion.bed.gz", sep="")
X <- fread(fn, stringsAsFactors=FALSE, header=T, data.table=F)
DP <- X[,5:ncol(X)]
rownames(DP) <- X$ID

##convariates
fn <- paste(prefix, "dispersionQTL/covariates/Tcell_PHA-DEX_covs_3PCs.txt", sep="")
cv.DP <- read.table(fn, header=TRUE, row.names=1)
names(cv.DP) <- gsub("\\.", "-", names(cv.DP))


## Gene and QTL
rn <- rownames(DP)
x <- grch38%>%filter(symbol=="RPS16")
ENSG <- rn[grepl(x[1,1], rn)]
varID <- dosages[grepl("19:39452746:TAA:T", dosages$varID), "varID"]    
symbol <- "RPS16"


### data.frame 1
y <- GE[ENSG,]
x <- dosages[varID, names(y)]
cv2 <- t(cv.mean[,names(y)])
d1 <- data.frame(y=as.numeric(y), x=round(as.numeric(x)),
   treats="PHA-DEX", index=1)
d1 <- cbind(d1, cv2)
reslm <- lm(y~0+factor(Batch)+Sex+as.numeric(cage1)+as.numeric(SCAIP1_6_genPC1)+
   as.numeric(SCAIP1_6_genPC2)+as.numeric(SCAIP1_6_genPC3)+
   as.numeric(PC1)+as.numeric(PC2)+as.numeric(PC3), data=d1)
d1$residuals <- as.numeric(reslm$residuals)
         
         
### data.frame 2
y <- DP[ENSG,]
x <- dosages[varID, names(y)]
cv2 <- t(cv.DP[,names(y)])
d2 <- data.frame(y=as.numeric(y), x=round(as.numeric(x)),treats="PHA-DEX", index=2)
d2 <- cbind(d2, cv2)
reslm <- lm(y~0+Batch+Sex+as.numeric(cage1)+as.numeric(SCAIP1_6_genPC1)+
   as.numeric(SCAIP1_6_genPC2)+as.numeric(SCAIP1_6_genPC3)+
   as.numeric(PC1)+as.numeric(PC2)+as.numeric(PC3), data=d2)
d2$residuals <- as.numeric(reslm$residuals)

###
## df <- rbind(d1[,c(1:4,19)], d2[,c(1:4,14)])
ref <- strsplit(varID, ":")[[1]][3]
alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])

p1 <- ggplot(data=d1,aes(x=factor(x), y=residuals))+
   geom_boxplot(
                fill="#1f78b4", outlier.shape=NA)+
   scale_x_discrete("19:39452746:TAA:T",
      labels=c("0"=paste(ref, "/", ref, sep=""),
      "1"=paste(ref, "/", alt, sep=""),
      "2"=paste(alt, "/", alt, sep="")))+
   geom_smooth(method='lm',  se=FALSE, colour="black")+
   geom_jitter(width=0.25, size=1)+
   facet_grid(.~index, labeller=labeller(index=c("1"="mean")))+ 
   ylab(bquote(~italic(.(symbol))~" normalized mean gene expression"))+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=10),
         strip.text.x=element_text(size=14))
###
## d2 <- d2[,c(1:4,14)]
p2 <- ggplot(data=d2,aes(x=factor(x), y=residuals))+
   geom_boxplot(
                fill="#1f78b4", outlier.shape=NA)+
   scale_x_discrete("19:39452746:TAA:T",
      labels=c("0"=paste(ref, "/", ref, sep=""),
      "1"=paste(ref, "/", alt, sep=""),
      "2"=paste(alt, "/", alt, sep="")))+
   geom_smooth(method='lm', se=FALSE, colour="black")+
   geom_jitter(width=0.25, size=1)+
   facet_grid(.~index, labeller=labeller(index=c("2"="variability")))+
   ylab(bquote(~italic(.(symbol))~" normalized gene expression variability"))+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_text(size=12),
         axis.text=element_text(size=10),
         strip.text.x=element_text(size=14))

figfn <- paste("./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.7_",
   symbol, "-", varID, ".png", sep="")
png(figfn, width=600, height=500, res=120)
print(plot_grid(p1, p2))
dev.off()


## df <- rbind(d1[,c("x","index","residuals")],
##             d2[,c("x", "index", "residuals")])

## p <- ggplot(data=df, aes(x=factor(x), y=residuals))+
##    geom_boxplot(fill="#1f78b4", outlier.shape=NA)+
##    scale_x_discrete("19:39452746:TAA:T",
##       labels=c("0"=paste(ref, "/", ref, sep=""),
##       "1"=paste(ref, "/", alt, sep=""),
##       "2"=paste(alt, "/", alt, sep="")))+
##    ## geom_smooth(data=df, aes(x=as.numeric(x), y=as.numeric(y)),
##    ##     method='lm', formula=y~as.numeric(x), se=FALSE)+
##    geom_jitter(width=0.25, size=1)+
##    facet_grid(.~index,
##               labeller=labeller(index=c("1"="mean", "2"="variability")),
##               scales="free")+ 
##    ylab(bquote(~italic(.(symbol))~" normalized expression mean (variability)"))+
##    ggtitle("T cell PHA+DEX")+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_text(size=12),
##          axis.text=element_text(size=10),
##          strip.text.x=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=12))

figfn <- paste("./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.7.2_",
   symbol, "-", varID, ".png", sep="")
png(figfn, width=550, height=500, res=120)
print(p)
dev.off()

   

########################
### LDA dynamic eQTL ###
########################

rm(list=ls())

option <- "DiagLDA2"
prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/"
col1 <- c("CTRL"="#828282", "LPS"="#fb9a99",
   "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###
fn <- paste(prefix, "dataset_contrast.txt", sep="")
dataset <- read.table(fn, header=F)
names(dataset) <- c("MCls", "LDA", "treat")

i <- 16
cell <- dataset[i,1]
lda <- dataset[i,2]
treat <- dataset[i,3]
col0 <- col1[[lda]]

### create directory
outdir <- "./6_DEG.CelltypeNew_output/Filter2_pub/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

###
## normalized GE files:
fn <- paste0(prefix, "LDA-eQTL_Julong/DiagLDA2/normalized_data/Tcell_ldaPHA-DEX_trtPHA-DEX.bed.gz", sep="")
data <- fread(fn, header=T, sep="\t", stringsAsFactors=FALSE, data.table=F)
GE <- data[,5:ncol(data)]
rownames(GE) <- data$ID

## txt file with dosages:
fn <- paste0(prefix,
   "LDA-eQTL_Julong/DiagLDA2/genotypes/dosages/1_eQTL_signif_dosages_Tcell_ldaPHA-DEX_trtPHA-DEX.txt.gz", sep="")
dosages <- fread(fn, sep="\t", stringsAsFactors=F, data.table=F)
rownames(dosages) <- dosages$varID

### LDA-eGene
fn <- paste(prefix,
   "LDA-eQTL_Julong/DiagLDA2/plots/eGene_SNP.16_Tcell.ldaPHA-DEX.trtPHA-DEX.txt", sep="")
res <- read.table(fn, header=TRUE)

## overlap with TWAS
## twas <- read.table("./DiagLDA2/lm_results/TWAS_immune-related-diseases.txt", header=TRUE)
## gene <- unique(twas$Gene)
## res2 <- res%>%filter(ENSG%in%gene)

###
res2 <- res%>%filter(symbol=="RCC1L")%>%
    arrange(interaction_ANOVA.pvalue)

## for (k in 1:nrow(res2)){
###
   k <- 1
   ENSG <- res2$ENSG[k]
   varID <- res2$varID[k]
   symbol <- res2$symbol[k]
   cat(ENSG, varID, "\n")
   y <- GE[ENSG,]
   x <- round(dosages[varID, names(y)])
   df <- data.frame(rn=names(y), y=as.numeric(y), x=as.numeric(x))%>%
      mutate(Bin=gsub(".*_", "", rn))
   ###
   ref <- strsplit(varID, ":")[[1]][3]
   alt <- gsub(";.*", "", strsplit(varID, ":")[[1]][4])
   snp <- gsub(".*;", "", varID)
   fig <- ggplot(data=df, aes(x=as.factor(x), y=y))+
      geom_boxplot(aes(alpha=Bin), fill=col0, outlier.shape=NA)+
  ## scale_fill_manual(values= c("0"="#4daf4a", "1"="#984ea3","2"="#ff7f00"))+
      scale_alpha_manual(values=c("1"=0.2, "2"=0.6, "3"=1), guide="none")+
      scale_x_discrete(snp, labels=c("0"=paste(ref, "/", ref, sep=""),
                                    "1"=paste(ref, "/", alt, sep=""),
                                    "2"=paste(alt, "/", alt, sep="")))+
      ggtitle("T cell PHA+DEX")+
      geom_jitter(width=0.25, size=1)+
      facet_grid(.~Bin, labeller=labeller(Bin=c("1"="1st tertile",
         "2"="2nd tertile","3"="3rd tertile"))) +
      ylab(bquote(italic(.(symbol))~ " normalized gene expresssion"))+
      theme_bw()+
      theme(axis.title.y=element_text(size=12),
            axis.text=element_text(size=10),
            plot.title=element_text(hjust=0.5,size=12),
            strip.text.x=element_text(size=12))
   ### 
   figfn <- paste(outdir, "/Figure1.8_", ENSG, "_", symbol, "_",
      varID, ".png", sep="")
   png(figfn, width=480, height=480, res=120)
   print(fig)
   dev.off()
   ###
   opfn <- paste(outdir, "1.8_LDA-eQTL.rds", sep="")
   write_rds(fig, opfn)

## }###



###
## p3 <- ggplot(df2)+
##    geom_density(aes(x, colour=factor(treat2)), size=1, fill=NA)+
##    scale_colour_manual(values=c("PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
##    xlab("LDA of PHA-DEX")+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.y=element_blank(),
##          axis.text=element_blank(),
##          axis.ticks=element_blank())
## figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure1.4.2_DLDA.png"
## png(figfn, width=300, height=250, res=120)
## print(p3)
## dev.off()

###
###
###
### ashg meeting

###
x <- seq(-4, 4, leng=100)
y <- dnorm(x, 0, 1)
d1 <- data.frame(x=x, y=y, treat="Control")
###
x <- seq(-2, 6,length=100)
y <- dnorm(x, 2, 1)
d2 <- data.frame(x=x, y=y,treat="Treatment")

df2 <- rbind(d1, d2)

p1 <- ggplot(df2)+
   geom_line(aes(x=x, y=y, colour=factor(treat)), size=1)+
   scale_color_manual(values=c("Control"="#4daf4a", "Treatment"="#ff7f00"))+
   ## ggtitle("Equal mean expression")+
   theme_void()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5,size=12)) 

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure3.1.0_eQTL.png"
png(figfn, width=400, height=250, res=120)
print(p1)
dev.off()


###
x <- rep(c(0,1,2),each=30)
y <- x*0.1+rnorm(90)
d1 <- data.frame(x=x, y=y, treat="Control")
##
y <- x*0.8+rnorm(90)*0.6
d2 <- data.frame(x=x, y=y, treat="Treatment")
df <- rbind(d1, d2)

p2 <- ggplot(data=df, aes(x=factor(x), y=y))+
   geom_boxplot(aes(fill=factor(x)), outlier.shape=NA)+
   scale_fill_manual("",
      values=c("0"="#4daf4a", "1"="#984ea3", "2"="#ff7f00"), guide="none")+
   scale_x_discrete(labels=c("0"="BB", "1"="AB","2"="AA"))+ 
   geom_jitter(width=0.15, size=0.5)+
   facet_grid(.~treat)+
   ylab("normalized gene expresssion")+
   theme_bw()+
   theme(legend.position="none",
         axis.title.y=element_text(size=12),
         axis.title.x=element_blank(),
         axis.text=element_text(size=10),
         strip.text.x=element_text(size=12))
figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure3.1.1_reQTL.png"
png(figfn, width=500, height=300, res=120)
print(p2)
dev.off()



###
###
x <- seq(-4, 4, leng=100)
y <- dnorm(x, 0, 1)
d1 <- data.frame(x=x, y=y, treat="Treatment")
###
x <- seq(-8,8,length=120)
y <- dnorm(x, 0, 2)
d2 <- data.frame(x=x, y=y,treat="Control")

df2 <- rbind(d1, d2)

p1 <- ggplot(df2)+
   geom_line(aes(x=x, y=y, colour=factor(treat)), size=1)+
   scale_color_manual(values=c("Control"="#4daf4a", "Treatment"="#ff7f00"))+
   ## ggtitle("Equal mean expression")+
   theme_void()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5,size=12)) 

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure3.2.0_vQTL.png"
png(figfn, width=350, height=250, res=120)
print(p1)
dev.off()

###
x <- rep(c(0, 1, 2),each=30)
y <- x*0.1+rnorm(90)+2
d1 <- data.frame(x=x, y=y, treat="Control")
##
y <- -x*0.5+rnorm(90)*0.5+3
d2 <- data.frame(x=x, y=y, treat="Treatment")
df <- rbind(d1, d2)

p2 <- ggplot(data=df, aes(x=factor(x), y=y))+
   geom_boxplot(aes(fill=factor(x)), outlier.shape=NA)+
   scale_fill_manual("",
      values=c("0"="#4daf4a", "1"="#984ea3", "2"="#ff7f00"), guide="none")+
   scale_x_discrete(labels=c("0"="BB", "1"="AB","2"="AA"))+ 
   geom_jitter(width=0.15, size=0.5)+
   facet_grid(.~treat)+
   ylab("normalized gene variability")+
   theme_bw()+
   theme(legend.position="none",
         axis.title.y=element_text(size=12),
         axis.title.x=element_blank(),
         axis.text=element_text(size=10),
         strip.text.x=element_text(size=12))
figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure3.2.1_rvQTL.png"
png(figfn, width=500, height=300, res=120)
print(p2)
dev.off()

###



