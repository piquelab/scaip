###
## library("rhdf5")
## library("corpcor")
## library(Matrix)
## library(MASS)
## library(scales)
library(tidyverse)
## library(parallel)
library(data.table)
## library(purrr)
## library(furrr)
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
library(ggrastr)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./9_RNA.dynamic2_output/Filter2_pub/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


#####################################
### Figure 4.1 and 4.2, LDA plots ###
#####################################
rm(list=ls())

df <- read_rds("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA2/1_MCls.Tcell.DLDA.rds")

###
col2 <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
###
###
df2 <- df%>%filter(treat2%in%c("CTRL", "LPS", "LPS-DEX"))%>%
   dplyr::select("z_LPS-DEX", "z_LPS", "treat2")%>%
   dplyr::rename("x"="z_LPS-DEX", "y"="z_LPS")

###
p1 <- ggplot(df2)+
   geom_point(aes(x=x, y=y, colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col2)+
   labs(x="DLDA of LPS+DEX", y="DLDA of LPS")+
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=10),
         axis.title=element_text(size=14))
p1 <- ggMarginal(p1, groupColour=TRUE, groupFill=FALSE, size=3)

## figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure6.1_DLDA.png"
## png(figfn, width=420, height=420, res=120)
## print(p1)
## dev.off()
##
## opfn <- "./9_RNA.dynamic2_output/Filter2_pub/6.1_DLDA.LPS.rds"
## write_rds(p1, opfn)

###
###
df2 <- df%>%filter(treat2%in%c("CTRL", "PHA", "PHA-DEX"))%>%
   dplyr::select("z_PHA-DEX", "z_PHA", "treat2")%>%
   dplyr::rename("x"="z_PHA-DEX", "y"="z_PHA")

###
p2 <- ggplot(df2)+
   geom_point(aes(x=x, y=y, colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col2)+
   labs(x="DLDA of PHA+DEX", y="DLDA of PHA")+
   theme_bw()+
   theme(legend.position="none",
         axis.text=element_text(size=10),
         axis.title=element_text(size=14))
p2 <- ggMarginal(p2, groupColour=TRUE, groupFill=FALSE, size=3)

## fig <- plot_grid(p1, p2, align="hv", axis="tb")
###
## figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure6.2_DLDA.png"
## png(figfn, width=420, height=420, res=120)
## print(p2)
## dev.off()


figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure6.1_comb.DLDA.png"
png(figfn, width=420, height=800, res=120)
plot_grid(p1, p2, nrow=2,
   align="hv", axis="lr")
   ## labels=NULL, label_x=0.05, label_y=0.9,
   ## label_size=18, label_fontface="plain")
dev.off()
###
## opfn <- "./9_RNA.dynamic2_output/Filter2_pub/6.2_DLDA.PHA.rds"
## write_rds(p2, opfn)



######################################################################
###  4.3 Dynamical changes of gene expression along the pseudotime ###
######################################################################

rm(list=ls())
### sliding window
slideFun <- function(X_sort, lda, win=0.1, step=0.001){
###
  win <- trunc(ncol(X_sort)*win)
  step <- trunc(ncol(X_sort)*step)
    
  X <- X_sort  
  nlen <- ncol(X) 

  Xnew <- NULL
  ldaNew <- NULL  
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    xi <- apply(X[,s0:s1], 1, sum)
    Xnew <- cbind(Xnew, xi)
    ldai <- mean(lda[s0:s1])
    ldaNew <- c(ldaNew, ldai)  
    s0 <- s0+step
  }
  res <- list(Xnew=Xnew, ldaNew=ldaNew)
  res  
}

###
getDynamicalMat <- function(X, meta, contrast, win=0.1, step=0.001){

 ##  treat0 <- contrast[1]
##   treat1 <- contrast[2]
## ###
##   barcode <- colnames(X)
##   x1 <- X[,grepl(treat0, barcode)]
##   mu1 <- apply(x1, 1, mean)
##   x2 <- X[,grepl(treat1, barcode)]
##   mu2 <- apply(x2, 1, mean)
##   diff <- data.frame(gene=gsub("S-|\\..*", "", rownames(X)), beta=mu2-mu1)%>%
##     arrange(beta) 
  ## geneSel <- diff$gene
  ## mat2 <- mat2[geneSel,]

### re-order cell
  meta_sort <- meta%>%arrange(z)

## treatment 1
  treat1 <- contrast[2]  
  meta1 <- meta_sort%>%filter(treats==treat1)  
  mat1 <- as.matrix(X[, meta1$NEW_BARCODE])
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp1 <- slideFun(mat1, lda=meta1$z, win=win, step=step)
  mat1 <- tmp1$Xnew
  lda1 <- tmp1$ldaNew  
### contrast treatment
  treat0 <- contrast[1]  
  meta0 <- meta_sort%>%filter(treats==treat0)  
  mat0 <- as.matrix(X[, meta0$NEW_BARCODE])
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  tmp0 <- slideFun(mat0, lda=meta0$z, win=win, step=step)
  mat0 <- tmp0$Xnew
  lda0 <- tmp0$ldaNew  
###
  mMax <- apply(cbind(mat1, mat0), 1, max)
    
  mat1 <- mat1[mMax>0,]
  mat0 <- mat0[mMax>0,]  
  mMax <- mMax[mMax>0]  
  mat1 <- sweep(mat1, 1, mMax, "/")
  mat1 <- as.matrix(mat1)
  ##
  mat0 <- sweep(mat0, 1, mMax, "/")
  mat0 <- as.matrix(mat0)
###
  mat_ls <- list(mat0=mat0, mat1=mat1, lda0=lda0, lda1=lda1)  
  mat_ls  
###    
}


###setting colors
col_fun <-  colorRamp2(seq(0, 1, len=11), rev(brewer.pal(11, "Spectral")))

res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")

##
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
dataset <- data.frame(MCls=rep(MCls, each=4), LDA=rep(LDA, times=4))

contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"),
   "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
coltreat <- c("CTRL"="#828282",
   "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")

option <- "DiagLDA2"
i <- 16
oneMCl <- dataset[i,1]
lda <- dataset[i,2]
contrast <- contrast_ls[[lda]]
treat0 <- contrast[1]
treat1 <- contrast[2] 
ii <- which(LDA==lda)
col0 <- coltreat[[treat1]]
###
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/3_MCl.",
   oneMCl, ".old.rds", sep="")
sc <- read_rds(fn)
X <- sc@assays$RNA@data
rn <- gsub("S-|\\..*", "", rownames(X))
rownames(X) <- rn
    
### sort cells by LDA
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/1_MCls.",
   oneMCl, ".DLDA.rds", sep="")
meta <- read_rds(fn)
meta2 <- meta%>%
   dplyr::select(NEW_BARCODE, treats, treat2)%>%mutate(z=meta[,ii])
res0 <- res%>%filter(MCls==oneMCl, contrast==lda, abs(beta)>0.5, qval<0.1)
    
###
### heatmap for treat1 by cluster genes
Xi <- X[as.character(res0$gene),]
mat_ls <- getDynamicalMat(Xi, meta2, contrast, win=0.1,step=0.001)
mat1 <- mat_ls[[2]]

### heatmap
### column annotation
z <- mat_ls$lda1
z2 <- (z-min(z))/(max(z)-min(z))
#z3 <- z2[z2>0.23] 
#breaks <- c(0,quantile(z3,probs=seq(0,1,length.out=99))) 
col2 <- colorRamp2(seq(0,1,length.out=100),
   colorRampPalette(c("white", coltreat[[treat1]]))(100))

col_ha <- HeatmapAnnotation(pseudotime=z2, col=list(pseudotime=col2),
   show_legend=FALSE,                         
   annotation_name_gp=gpar(fontsize=8),                         
   simple_anno_size=unit(0.3,"cm") )
###
# <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
b <- res0$beta
th0 <- as.numeric(quantile(abs(b),probs=0.99)) 
b2 <- b[abs(b)<th0] ### quantile(abs(b),probs=0.99)
breaks <- c(min(b),quantile(b2,probs=seq(0,1,length.out=98)),max(b)) 
col2 <- colorRamp2(breaks,
   colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) ) 
row_ha <- rowAnnotation(LFC=res0$beta, col=list(LFC=col2),
   annotation_legend_param=list(LFC=list(grid_width=unit(0.3,"cm"),
   labels_gp=gpar(fontsize=8), title_gp=gpar(fontsize=10), title="LFC")),    
   width=unit(0.5,"cm"), annotation_name_gp=gpar(fontsize=8)) 

###
fig1 <- Heatmap(mat1, col=col_fun,
   cluster_rows=TRUE, cluster_columns=FALSE,
   row_km=4, show_parent_dend_line=FALSE,
   show_row_names=FALSE, show_column_names=FALSE,
   column_title=paste("pseudotime(", treat1, ")", sep=""),
   column_title_gp=gpar(fontsize=10),
   top_annotation=col_ha, right_annotation=row_ha,
   heatmap_legend_param=list(title="Expression",
   title_gp=gpar(fontsize=10),
   labels_gp=gpar(fontsize=8), grid_width=unit(0.3, "cm")),
   use_raster=TRUE, raster_device="png")


### output directory
outdir <- paste("./9_RNA.dynamic2_output/Filter2_pub/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T) 

figfn <- paste(outdir, "Figure6.3_16_Tcell_ldaPHA-DEX_trtPHA-DEX.png", sep="")
png(figfn, height=800, width=650, res=150)
set.seed(0)
fig1 <- draw(fig1)
r.list <- row_order(fig1)
r.dend <- row_dend(fig1) 
dev.off()

opfn <- paste(outdir, "6.3_LDA.dynamic.rds", sep="")
write_rds(fig1, opfn)


## fn <- paste(outdir, "6.3_LDA.dynamic.rds", sep="")
## fig <- read_rds(fn)
## figfn <- paste(outdir, "Figure6.3_16_Tcell_ldaPHA-DEX_trtPHA-DEX.png", sep="")
## png(figfn, height=600, width=520, res=120)
## set.seed(0)
## print(fig)# <- draw(fig)
## dev.off()


###cluster genes
clu_df <- lapply(names(r.list),function(i){
   out <- data.frame(genes=rownames(mat1)[r.list[[i]]],
                     Cluster=i, stringsAsFactors=FALSE)
   out
})%>%do.call(rbind,.)
rownames(clu_df) <- clu_df$gene
###
opfn <- paste(outdir, "6.3_Tcell_ldaPHA-DEX_trtPHA-DEX.cluster.csv", sep="")
write.csv(clu_df, opfn, row.names=FALSE)


############################################
### figure 4.4, average gene expression  ###
############################################
###
Xi <- X[as.character(res0$gene),]
mat_ls <- getDynamicalMat(Xi, meta2, contrast, win=0.1, step=0.001)
mat1 <- mat_ls[[2]]
newLDA <- mat_ls$lda1
col0 <- coltreat[["PHA-DEX"]]

fn <- paste(outdir, "6.3_Tcell_ldaPHA-DEX_trtPHA-DEX.cluster.csv", sep="")
clu_gene <- read.csv(fn)

df2 <- map_dfr(1:4, function(i){
   clu_gene2 <- clu_gene%>%filter(Cluster==i)  
   gene0 <- clu_gene2$genes#[1:10]     
   y <- colMeans(mat1[gene0,])
   ii <- 5-i
   cvt <- data.frame(z=newLDA, y=y, cluster=i, cluster2=letters[ii])
   cvt
})

###
### (2) plot curve
df2 <- df2%>%drop_na(y)
fig <-  ggplot(df2, aes(x=z, y=y))+
   geom_smooth(method="loess", colour=col0, se=F)+
   ylab("Relative Expression")+
   xlab("pseudotime")+
   ylim(c(0,1))+
   facet_wrap(vars(cluster2), nrow=4,
      labeller=labeller(cluster2=c("a"="cluster 4",
      "b"="cluster 3", "c"="cluster 2", "d"="cluster 1")))+ 
   theme_bw()+
   theme(#plot.title=element_text(hjust=0.5, size=10),
         axis.title.x=element_text(size=12),
         axis.title.y=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_text(size=10),
         axis.ticks.x=element_blank(),
         strip.text=element_text(size=12))

###
## fig2 <- plot_grid(fig_ls[[4]], fig_ls[[3]], fig_ls[[2]], fig_ls[[1]], nrow=4)
figfn <- paste(outdir, "Figure6.4_fitting.png", sep="")
png(figfn, width=280, height=800, res=120)
print(fig)
dev.off()

##
opfn <- paste(outdir, "6.4_fitting.rds", sep="")
write_rds(fig, opfn)

## fn <- paste(outdir, "6.4_fitting.rds", sep="")
## fig <- read_rds(fn)+theme(axis.text.y=element_text(size=8))
## figfn <- paste(outdir, "Figure6.4_average.fitting.png", sep="")
## png(figfn, width=250, height=600, res=120)
## print(fig)
## dev.off()


########################################
### Figure 7, example dynamic eQTL ###
########################################


################
### boxplots ###
################

dataset <- read.table("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/LDA-eQTL_Julong/dataset_contrast.txt",header=F)

col1 <- c("CTRL"="#828282",
   "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

i <- 16
cell <- dataset[i,1]
lda <- dataset[i,2]
treat <- dataset[i,3]
col0 <- col1[[lda]]

###
## normalized GE files:
prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/LDA-eQTL_Julong/DiagLDA2/"
fn <- paste0(prefix, "1_normalized.data/", cell, "_lda", lda, "_trt", treat, ".bed.gz")
data <- fread(fn, header=T, sep="\t", stringsAsFactors=FALSE, data.table=F)
GE <- data[,5:ncol(data)]
rownames(GE) <- data$ID

## txt file with dosages:
fn <- paste0(prefix, "3_genotypes/dosages/1_eQTL_signif_dosages_", cell, "_lda", lda, "_trt", treat, ".txt.gz")
dosages <- fread(fn, sep="\t", stringsAsFactors=F, data.table=F)
rownames(dosages) <- dosages$varID

### LDA-eGene
fn <- paste(prefix, "examples/eGene_SNP.", i, "_", cell, ".lda",
   lda, ".trt", treat, ".txt", sep="")
res <- read.table(fn, header=TRUE)

## overlap with TWAS

### candidate gene list
## df.gene <- data.frame(symbol=c("CSF3", "ABO", "HLA-DRB5", "HLA-DRB5", "HLA-DRB5",  "GADD45G",
##                                "PHRF1", "RGCC", "TRAF1", "TRAF1"),
##    rs=c("rs35938904", "rs34426871", "rs71536582", "rs116611418", "6:32487309:T:C",  "rs11265809",
##         "rs11246184", "rs367833719", "rs10985112", "rs60692488"))

df.gene <- data.frame(symbol=c("ABO", "HLA-DRB5",  "GADD45G"),
   rs=c("rs34426871", "rs116611418",  "rs11265809"))

##
## df.gene <- df.gene%>%dplyr::filter(symbol=="HLA-DRB5")
boxs <- lapply(1:nrow(df.gene), function(ii){

   symbol0 <- df.gene[ii,1]
   rs0 <- df.gene[ii,2]
   res2 <- res%>%filter(symbol==symbol0, grepl(rs0, varID)) ### HLRC1
   k <- 1
   ENSG <- res2$ENSG[k]
   varID <- res2$varID[k]
   symbol <- res2$symbol[k]  

y <- GE[ENSG,]
x <- round(dosages[varID, names(y)])
df <- data.frame(rn=names(y), y=as.numeric(y), x=as.numeric(x))%>%
   mutate(Bin=gsub(".*_", "", rn))    

ref <- strsplit(varID, ":")[[1]][3]
alt <- gsub(";.*", "", strsplit(varID, ":")[[1]][4])
rs <- gsub(".*;", "", varID)  
p0 <- ggplot(data=df, aes(x=as.factor(x), y=y))+
   geom_boxplot(aes(alpha=Bin), fill=col0, outlier.shape=NA)+
   ## scale_fill_manual(values= c("0"="#4daf4a", "1"="#984ea3","2"="#ff7f00"))+
   scale_alpha_manual(values=c("1"=0.2, "2"=0.6, "3"=1), guide="none")+
   scale_x_discrete(rs, labels=c("0"=paste(ref, "/", ref, sep=""),
      "1"=paste(ref, "/", alt, sep=""),
      "2"=paste(alt, "/", alt, sep="")))+
   geom_smooth(method='lm', se=F, color="black", size=0.6, aes(group=Bin))+
   geom_jitter(width=0.25, size=0.8)+
   facet_grid(.~Bin, labeller=labeller(Bin=c("1"="1st tertile",
      "2"="2nd tertile", "3"="3rd tertile"))) +
   ## ylab(bquote(italic(.(symbol))~ " normalized gene expresssion"))+
   ylab("Normalized gene expression")+ 
   ggtitle(bquote("T cell PHA+DEX"~"("~italic(.(symbol))~")"))+
   theme_bw()+
   theme(axis.title.x=element_text(size=10),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(size=9),
         axis.text.y=element_text(size=10),
         strip.text.x=element_text(size=10),
         plot.title=element_text(hjust=0.5,size=12))
       
p0
})    
    

## outdir <- "./9_RNA.dynamic2_output/Filter2_pub/"
## figfn <- paste(outdir, "/Figure6.7_", ENSG, ".", symbol, ".", varID, ".png", sep="")
## png(figfn, width=450, height=550, res=120)
## print(fig)
## dev.off()

## ###
## outdir <- "./9_RNA.dynamic2_output/Filter2_pub/"
## opfn <- paste(outdir, "/6.7_LDA-eQTL.rds", sep="")
## write_rds(fig, opfn)


#######################
### scatter plots   ###
#######################


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
    ## tmp <- d%>%filter(dosage==2)  
    ## cat(s0, s1, tmp$y, "\n")  
    s0 <- s0+step      
  }
  #cvt <- cbind(cvt, cvt2)
  cvt2
}

###
### counts data
sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
X <- sc@assays$RNA@data
X <- X[grepl("S-",rownames(X)),]
rn <- gsub("S-|\\..*", "", rownames(X))
rownames(X) <- rn
meta <- sc@meta.data%>%dplyr::select(NEW_BARCODE, BEST.GUESS)


###
### read dosage
prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/LDA-eQTL_Julong/"
genfn <- paste(prefix, "check_output/SCAIP1-6_filtered.vcf.txt.gz", sep="")
dosages <- fread(genfn, data.table=F, stringsAsFactors=F)
fn <- paste(prefix, "check_output/sample.txt", sep="")
sample <- read.table(fn)
names(dosages) <- c("chr", "pos", "varID", sample$V1)
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID


###
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
treats <- c("LPS-EtOH", "LPS-DEX", "PHA-EtOH", "PHA-DEX")
dataset <- data.frame(MCls=rep(MCls, each=4),
   LDA=rep(LDA, times=4), treats=rep(treats, times=4))

i <- 16
cell <- dataset[i,1]
lda <- dataset[i,2]
treat <- dataset[i,3]

### read LDA-interaction eQTL
fn <- paste(prefix, "DiagLDA2/examples/eGene_SNP.", i, "_", cell, ".lda", lda, ".trt", treat, ".txt", sep="")
signif <- read.table(fn, header=T)
## pick up SNP-gene pairs

###
###
## df.gene <- data.frame(symbol=c("KLRC1", "HLA-DQA1", "IFNG"),
##                       rs=c("rs140718479", "rs9271790", "rs1362463"))

## df.gene <- data.frame(symbol=c("CSF3", "ABO", "HLA-DRB5", "HLA-DRB5", "HLA-DRB5",  "GADD45G",
##                                "PHRF1", "RGCC", "TRAF1", "TRAF1"),
##    rs=c("rs35938904", "rs34426871", "rs71536582", "rs116611418", "6:32487309:T:C",  "rs11265809",
##         "rs11246184", "rs367833719", "rs10985112", "rs60692488"))

##

df.gene <- data.frame(symbol=c("ABO", "HLA-DRB5",  "GADD45G"),
   rs=c("rs34426871", "rs116611418",  "rs11265809"))

points <- lapply(1:nrow(df.gene), function(ii){

   symbol0 <- df.gene[ii,1]
   rs0 <- df.gene[ii,2]
   signif2 <- signif%>%filter(symbol==symbol0, grepl(rs0, varID)) ### HLRC1
## signif2 <-signif%>%filter(symbol=="HLA-DQA1") ### KLRC1

k <- 1
ENSG <- signif2$ENSG[k] 
ENSG2 <- signif2$ENSG2[k]
symbol <- signif2$symbol[k]
varID <- signif2$varID[k]
gen <- dosages[varID,-(1:3)]
gen <- round(as.numeric(gen[1,]))
names(gen) <- names(dosages)[-(1:3)]

###
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA2/",
   "/1_MCls.", cell, ".DLDA.rds", sep="")
meta2 <- read_rds(fn)%>%left_join(meta, by="NEW_BARCODE")
cvt <- meta2%>%dplyr::select(NEW_BARCODE, treats, BEST.GUESS)%>%
   mutate(z=meta2[,4], y=X[ENSG2,NEW_BARCODE], dosage=gen[meta2$BEST.GUESS])

###
## cvt2 <- map_dfr(c("PHA-EtOH", "PHA-DEX"), function(oneX){
##    ###
##    cvt0 <- cvt%>%filter(treats=="PHA-DEX")
##    cvt0 <- slideFun2(cvt0, win=0.2, step=0.01)
##    cvt0 <- as.data.frame(cvt0)%>%mutate(treats=oneX)
##    cvt0
## })

cvt0 <- cvt%>%dplyr::filter(treats=="PHA-DEX")
cvt2 <- slideFun2(cvt0, win=0.2, step=0.01)    
cvt2 <- cvt2%>%drop_na(y)
cvt2$y2 <- cvt2$y/max(cvt2$y)
    
### (2) plot curve
ref <- strsplit(varID, ":")[[1]][3]
alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])

p0 <- ggplot(cvt2, aes(x=z_ave, y=y2))+
   geom_point(aes(colour=factor(dosage)), size=0.8)+  
   geom_smooth(aes(colour=factor(dosage)), method="loess",
       size=0.5, se=F)+
   scale_colour_manual("genotype",
        values=c("0"="#7570b3", "1"="#1b9e77", "2"="#d95f02"),
        ## values=c("0"="#4daf4a", "1"="#984ea3", "2"="#ff7f00"),               
        ## labels=c("0"="BB","1"="AB","2"="AA"))+
        labels=c("0"=paste(ref, "/", ref, sep=""),
                 "1"=paste(ref, "/", alt, sep=""),
                 "2"=paste(alt, "/", alt, sep="")))+ 
    
   scale_y_continuous(limits=c(0,max(cvt2$y2)))+
   xlab("LDA of PHA+DEX")+
   ## xlab("LDA pseudotime")+ 
   ylab("Relative changes")+ 
   ## ylab(bquote(~"Relative changes ("~italic(.(symbol))~")"))+      
   ## ggtitle("Tcell PHA+DEX")+
   theme_bw()+
   theme(## plot.title=element_text(hjust=0.5, size=12),
         ## legend.position=c(0.3,0.85),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.key=element_blank(),
         legend.key.size=grid::unit(1.2,"lines"),
         axis.title.x=element_text(size=10),
         axis.title.y=element_text(size=12),
         axis.text=element_text(size=10))
##
if (ii==1){
    p0 <- p0+theme(legend.position=c(0.3,0.85))
}else{
    p0 <- p0+theme(legend.position=c(0.7,0.85))  
}    
p0
})    

###
## fig1 <- plot_grid(boxs[[1]], points[[1]],
##    nrow=2, ncol=1,
##    align="v", axis="lr")
##    ## label_x=0.05,
##    ## labels=c("A", ""), label_fontface="plain")

## fig2 <- plot_grid(boxs[[2]], points[[2]],
##    nrow=2, ncol=1,
##    align="v", axis="lr")## ,
##    ## label_x=0.05,
##    ## labels=c("B", ""), label_fontface="plain")

## fig3 <- plot_grid(boxs[[3]], points[[3]],
##    nrow=2, ncol=1,
##    align="v", axis="lr")## ,
   ## label_x=0.05,
   ## labels=c("C", ""), label_fontface="plain")

## figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure7.3_comb.png"
## png(figfn, width=450, height=680, res=120)
## fig3
## dev.off()

## figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure7.0_comb.pdf"
## pdf(figfn, width=10, height=7)
## plot_grid(fig1, fig2, fig3,
##    nrow=1, ncol=3,       
##    align="h", axis="bt", labels=NULL)
## dev.off()

## ##
## figfn <- "./9_RNA.dynamic2_output/Filter2_pub/Figure8.1_KLRC1.fitting.png"
## png(figfn, width=480, height=480, res=120)
## points[[1]]
## dev.off()

### 2
symbol <- df.gene$symbol
for (i in 1:length(symbol)){
   symbol0 <- symbol[i]
   cat(i, symbol0, "\n")
   ###
   fig0 <- plot_grid(boxs[[i]],points[[i]], nrow=2, ncol=1, align="v", axis="lr")   
   figfn <- paste("./9_RNA.dynamic2_output/Filter2_pub/Figure7",  "_", symbol0, ".", i, ".png", sep="")
   png(figfn, width=400, height=820, res=120)
   print(fig0)
   dev.off()
}
