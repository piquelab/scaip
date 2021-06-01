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
library(SeuratDisk)
library(harmony)
library(annotables) 
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

outdir <- "./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA3/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


############################################################
### 2021-05-03, psedotime definition based on 6,571 DEGs ###
###       diagonal linear Discriminant analysis          ###
###       Last modified by Julong wei, 2021-05-03        ###
############################################################

rm(list=ls())

########################
### 1. calculate LDA ### 
########################

#res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/3_Batch1456.meta.rds")
res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
res <- res%>%mutate(gamma=ifelse( (abs(beta)>0.5)&(qval<0.1), 1, 0))
##
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
option <- "DiagLDA3"
tmp <- mclapply(1:4, function(i){
###
   oneMCl <- MCls[i]
   s1 <- Sys.time()
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/3_MCl.", oneMCl, ".old.rds", sep="")
   sc <- read_rds(fn)

   X <- sc@assays$RNA@data
   rn <- gsub("S-|\\..*", "", rownames(X))
   rownames(X) <- rn 

   meta <- sc@meta.data
   treats <- unique(meta$treats)
   ncell <- ncol(X)
   ntreat <- length(treats)
    
   ### variance across genes  
   ## Xc <- map_dfc(treats, function(one){
   ##    ii <- as.character(meta[meta$treats==one, "NEW_BARCODE"])
   ##    Xi <- X[,ii]
   ##    mu_k <- apply(Xi, 1, mean)
   ##    xc <- sweep(Xi, 1, mu_k, "-")     
   ##    rowSums(xc*xc)      
   ## })
   mu <- apply(X, 1, mean)
   Xc <- sweep(X, 1, mu, "-") 
   var <- rowSums(Xc*Xc)/(ncell-1)
   names(var) <- rn 

   ###
   contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                       "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                       "PHA"=c("CTRL", "PHA-EtOH"),
                       "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
   contrast_nn <- names(contrast_ls)
   metaNew <- sapply(contrast_nn, function(nn){
      oneX <- contrast_ls[[nn]] 
      ### extract lfc as weights
      res0 <- res%>%filter(MCls==oneMCl, contrast==nn)
      gamma <- rep(0, length(rn))
      names(gamma) <- rn
      gamma[res0$gene] <- 1

      ## res0 <- res%>%filter(MCls==oneMCl, contrast==nn)
      ## lfc <- abs(res0$beta)
      ## lfc[is.na(lfc)] <- 0
      ## names(lfc) <- as.character(res0$gene)
      ## lfc0 <- lfc[rn]
      
      ### matrix 1 
      ii <- meta%>%filter(treats==oneX[1])%>%dplyr::pull(NEW_BARCODE)
      X1 <-X[,ii]
      mu1 <- apply(X1, 1, mean)

      ### matrix 2
      ii <- meta%>%filter(treats==oneX[2])%>%dplyr::pull(NEW_BARCODE)
      X2 <-X[,ii]
      mu2 <- apply(X2, 1, mean)

      Diff <- as.matrix((mu2-mu1)*(1/var)) ##

      if ( option=="DiagLDA"){
         s0 <- rep(1, length(rn))
      }    
      ##
      if( option=="DiagLDA2"){
         gamma <- res0$gamma
         names(gamma) <- res0$gene
         gamma[is.na(gamma)] <- 0
         s0 <- gamma[rn]
      }
      ###
      if (option=="DiagLDA3"){
         lfc <- abs(res0$beta)
         names(lfc) <- res0$gene
         lfc[is.na(lfc)] <- 0
         s0 <- lfc[rn]
      }    

      Diff <- as.matrix((mu2-mu1)*(1/var)*s0)
      Diff[is.na(Diff)] <- 0
      ##
      ## X12 <- cbind(X1, X2)
      ## mu12 <- apply(X12, 1, mean)
      ## X12 <- sweep(X12, 1, mu12, "-")
      z1 <- as.vector(crossprod(Xc, Diff))
      z1
   })
   metaNew <- as.data.frame(metaNew) 
   names(metaNew) <- paste("z_", contrast_nn, sep="")  
   metaNew <- metaNew%>%
      mutate(NEW_BARCODE=colnames(X),
             treats=gsub("^SCAIP.{1,3}-|_.*", "", NEW_BARCODE),
             treat2=gsub("-EtOH", "", treats)) 
 
   ### 
   opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/",
                 option, "/1_MCls.", oneMCl, ".DLDA.rds", sep="")
   write_rds(metaNew, opfn)
      
   s2 <- Sys.time()
   d12 <- difftime(s2, s1, units="mins")
   cat(oneMCl, ":", d12, "\n")
}, mc.cores=1)



#########################
### summary LDA plots ###
#########################
### scatter plots plus density plots

LDAplot <- function(df, labXY=c("x","y") ){
    
  coltreat <- c("CTRL"="#828282",
                "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
    
  MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

  fig_ls <- lapply(MCls, function(oneMCl){
    df2 <- df%>%filter(MCls==oneMCl)
    fig0 <- ggplot(df2, aes(x=x, y=y))+
      geom_point(aes(colour=factor(treat2)), size=0.1)+
      scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=2    )))+
    labs(x=labXY[1], y=labXY[2])+
    ggtitle(oneMCl)+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=10),
          axis.title=element_text(size=6),
          axis.text=element_text(size=6))
    fig0 <- ggMarginal(fig0, groupColour=T, groupFill=F, size=2)
    fig0
  })

  lab2 <- c("CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS+DEX",
            "PHA"="PHA", "PHA-DEX"="PHA+DEX")
  legend2 <- get_legend(
    ggplot(df%>%filter(MCls=="Bcell"), aes(x, y))+
    geom_point(aes(colour=factor(treat2)), size=0.1)+
    scale_colour_manual(values=coltreat, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.background=element_rect(colour=NA, fill=NA),
          legend.text=element_text(size=8),
          legend.key=element_rect(fill=NA),
          legend.key.size=grid::unit(1,"lines"))
   ) ##end for legend2

  figures <- list(fig_ls=fig_ls, fig_legend=legend2)

  figures

} ###    
    
    

### Read data
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
option <- "DiagLDA3"
df <- map_dfr(MCls, function(oneMCl){
   ##
   cat(oneMCl, "\n")
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
               "/1_MCls.", oneMCl, ".DLDA.rds", sep="")
   dd <- read_rds(fn)%>%mutate(MCls=oneMCl)
   dd
})

### LPS vs PHA
df2 <- df%>%
   dplyr::select(NEW_BARCODE,treats, treat2, MCls, z_LPS, z_PHA)%>%
   dplyr::rename(x=z_LPS, y=z_PHA)
figures <- LDAplot(df2,labXY=c("LDA_LPS", "LDA_PHA"))
fig_ls <- figures$fig_ls
legend2 <- figures$fig_legend

figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure1.X1_LPSandPHA.png", sep="")
png(figfn, width=900, height=700, res=130)
fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


### LPS-DEX vs PHA-DEX
df2 <- df%>%
   dplyr::select(NEW_BARCODE, treats, treat2, MCls, "z_LPS-DEX", "z_PHA-DEX")%>%
   dplyr::rename(x="z_LPS-DEX", y="z_PHA-DEX")
figures <- LDAplot(df2, labXY=c("LDA_LPS-DEX", "LDA_PHA-DEX"))
fig_ls <- figures$fig_ls
legend2 <- figures$fig_legend

figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure1.X2_LPS-DEXandPHA-DEX.png", sep="")
png(figfn, width=900, height=700, res=130)
fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


### LPS vs LPS-DEX
df2 <- df%>%
   dplyr::select(NEW_BARCODE, treats, treat2, MCls, z_LPS, "z_LPS-DEX")%>%
   dplyr::rename(x=z_LPS, y="z_LPS-DEX")
figures <- LDAplot(df2, labXY=c("LDA_LPS", "LDA_LPS-DEX"))
fig_ls <- figures$fig_ls
legend2 <- figures$fig_legend

figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure1.X3_LPSandLPS-DEX.png", sep="") 
png(figfn, width=900, height=700, res=130)
fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


### PHA vs PHA-DEX
df2 <- df%>%
   dplyr::select(NEW_BARCODE, treats, treat2, MCls, z_PHA, "z_PHA-DEX")%>%
   dplyr::rename(x=z_PHA, y="z_PHA-DEX")
figures <- LDAplot(df2, labXY=c("LDA_PHA", "LDA_PHA-DEX"))
fig_ls <- figures$fig_ls
legend2 <- figures$fig_legend
figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure1.X4_PHAandPHA-DEX.png", sep="")
png(figfn, width=900, height=700, res=130)
fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


#######################################
### split Tcells into CD4+ and CD8+ ###
#######################################

## sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
## x <- sc@meta.data%>%dplyr::select(NEW_BARCODE, seurat_clusters)
## fn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/1_MCls.Tcell.DLDA.rds"
## dd <- read_rds(fn)%>%mutate(MCls="Tcell")%>%
##    left_join(x, by="NEW_BARCODE")

## dd$MCl2 <- "CD4+"
## dd$MCl2[dd$seurat_clusters==5] <- "CD8+"

## df2 <- dd
## MCl2 <- c("CD4+", "CD8+")
## fig_ls <- lapply(MCl2, function(oneMCl){
##    df0 <- df2%>%filter(MCl2==oneMCl)
##    fig0 <- ggplot(df0, aes(x=z_LPS, y=z_PHA))+
##       geom_point(aes(colour=factor(treat2)), size=0.1)+
##       scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
##       ggtitle(oneMCl)+
##       theme_bw()+
##       theme(legend.position="none",
##             plot.title=element_text(hjust=0.5, size=10),
##             axis.title=element_text(size=6),
##             axis.text=element_text(size=6))
##    fig0 <- ggMarginal(fig0, groupColour=T, groupFill=F, size=2)
##    fig0
## })

## lab2 <- c("CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS+DEX",
##                      "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## legend2 <- get_legend(
##    ggplot(df2%>%filter(MCl2=="CD8+"), aes(z_LPS_DEX, z_PHA_DEX))+
##    geom_point(aes(colour=factor(treat2)), size=0.1)+
##    scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.background=element_rect(colour=NA, fill=NA),
##          legend.text=element_text(size=8),
##          legend.key=element_rect(fill=NA),
##          legend.key.size=grid::unit(0.8,"lines"))
##    )

## ##
## png("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/Figure1.X1.Tcell.png", width=800, height=400, res=130)
## fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]], nrow=1, ncol=2, align="hv",axis="tb")
## print(plot_grid(fig1, legend2, rel_widths=c(5,1)))
## dev.off()
## ###
## ###


################################################
### 2. defined Three bins according to 4 LDAs ###
################################################

rm(list=ls())

option <- "DiagLDA2"

outdir <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/DLDA_Bin/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=TRUE)

###Fun, class Bin   
Binfun <- function(LDA){
   x <- quantile(LDA,probs=c(1/3, 2/3))
   Bin <- rep(1,length(LDA))
   Bin[(LDA>=x[1]&LDA<x[2])] <- 2 
   Bin[LDA>=x[2]] <- 3
   Bin
}

### class bin
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){
###  
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/1_MCls.", oneMCl, ".DLDA.rds", sep="")
   meta <- read_rds(fn)
   treats <- sort(unique(meta$treats)) 

   meta <- map_dfr(treats, function(ii){
      meta0 <- meta%>%filter(treats==ii)
      meta0$Bin_LPS <- Binfun(meta0$z_LPS)
      meta0$"Bin_LPS-DEX" <- Binfun(meta0$"z_LPS-DEX")
      meta0$Bin_PHA <- Binfun(meta0$z_PHA)
      meta0$"Bin_PHA-DEX" <- Binfun(meta0$"z_PHA-DEX")
      meta0
   })
###
   opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
                 "/DLDA_Bin/1_MCls.", oneMCl, ".Bin.rds", sep="")
   write_rds(meta, opfn)
   cat(oneMCl, 0, "\n")
}


##########################################
### 3. Average GE for each combination ###
##########################################

rm(list=ls())

############################
### (1). get counts data ###
############################

sc <- read_rds("./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds")
count <- sc@assays$RNA@counts         
x <- sc@meta.data%>%dplyr::select(NEW_BARCODE, BEST.GUESS, MCls, BATCH)

grch38_unq <- grch38%>%
              distinct(ensgene, .keep_all=T)%>%
              dplyr::select(ensgene, symbol, chr, biotype)
vars <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn),
               ensgene2=gsub("[SU]-", "", rn), 
               uns=grepl("S-",rn), rnz=rowSums(count))%>%
        left_join(grch38_unq, by="ensgene")

autosome <- as.character(1:22)        
varsSel <- vars%>%filter(uns, rnz>20, chr%in%autosome, grepl("protein_coding", biotype))          

count <- count[varsSel$rn,]
rownames(count) <- varsSel$ensgene2

#####################################
### (2), average expression value ###
#####################################

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Bins <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
option <- "DiagLDA2"
###
for (oneMCl in MCls){

fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/DLDA_Bin/1_MCls.", oneMCl, ".Bin.rds", sep="")
meta0 <- read_rds(fn)
x0 <- x%>%filter(NEW_BARCODE%in%meta0$NEW_BARCODE)
meta0 <- meta0%>%left_join(x0, by="NEW_BARCODE")
count0 <- count[, meta0$NEW_BARCODE]

### Bin    
for(i in 1:4){
   ii <- i+7    
   bti <- meta0%>%transmute(bti=paste(BEST.GUESS, MCls, treats, BATCH, meta0[,ii], sep="_")) %>% unlist %>% factor
   ncell <- data.frame(x=bti)%>%
         group_by(x)%>%
         summarise(ncell=n(),.groups="drop")
   X <- model.matrix(~0+bti)
   colnames(X) <- gsub("bti","",colnames(X))
   rownames(X) <- meta0$NEW_BARCODE

   x1 <- as.character(ncell$x)
   x2 <- as.character(colnames(X))
   x3 <- as.character(colnames(count0))
   x4 <- as.character(rownames(X))
   cat(oneMCl, "bti", identical(x1,x2), "\n")
   cat(oneMCl, "cells", identical(x3,x4), "\n")
#
####
   YtX <- count0 %*% X 
   YtX <- as.matrix(YtX)

###
   prefix <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
                   "/DLDA_Bin/2_", oneMCl, ".", Bins[i], sep="")
   opfn1 <- paste(prefix, "_YtX.sum.rds", sep="")
   write_rds(YtX, opfn1)

###    
   YtX_ave <- sweep(YtX, 2, ncell$ncell, "/")
   opfn2 <- paste(prefix, "_YtX.ave.rds", sep="")
   write_rds(YtX_ave, opfn2)

###    
   opfn3 <- paste(prefix, "_ncell.rds", sep="")
   write_rds(ncell, opfn3)
    
} ### Bin

} ### MCls


#################################################################
### Dynamical changes of gene expression along the pseudotime ###
#################################################################

rm(list=ls())
### sliding window
slideFun <- function(X_sort, win=0.2, step=0.005){
###
  win <- trunc(ncol(X_sort)*win)
  step <- trunc(ncol(X_sort)*step)
    
  X <- X_sort  
  nlen <- ncol(X) 

  Xnew <- NULL
  s0 <- 1
  while(TRUE){
    ##
    s1 <- s0+win-1
    if (s1>nlen) break
    xi <- apply(X[,s0:s1], 1, sum)
    Xnew <- cbind(Xnew, xi)
    s0 <- s0+step
  }
  Xnew  
}

###
getDynamicalMat <- function(X, meta, contrast, win=0.2, step=0.005){

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
  mat1 <- slideFun(mat1, win=win, step=step)
### contrast treatment
  treat0 <- contrast[1]  
  meta0 <- meta_sort%>%filter(treats==treat0)  
  mat0 <- as.matrix(X[, meta0$NEW_BARCODE])
  #rownames(mat2) <- gsub("S-|\\..*", "", rownames(mat2))
  mat0 <- slideFun(mat0, win=win, step=step)
    
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
  mat_ls <- list(mat0, mat1)  
  mat_ls  
###    
}


###setting colors
col_fun <-  colorRamp2(seq(0, 1, len=11), rev(brewer.pal(11, "Spectral")))

res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/2_meta.rds")

##
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
contrast_LDA <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"),
   "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))

option <- "DiagLDA2"
### loop by cell type
for (oneMCl in MCls){
    
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

###
   LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
   ### loop by LDA 
   for (i in 1:4){
      lda <- LDA[i]
      contrast <- contrast_LDA[[lda]]
      meta2 <- meta%>%
         dplyr::select(NEW_BARCODE, treats, treat2)%>%mutate(z=meta[,i])

      res0 <- res%>%filter(MCls==oneMCl, contrast==lda, abs(beta)>0.5, qval<0.1)
    
###
### heatmap for treat1 by cluster genes
      treat1 <- contrast[2]
      ## meta3 <- meta2%>%filter(treats==treat1)
      Xi <- X[as.character(res0$gene),]
      mat_ls <- getDynamicalMat(Xi, meta2, contrast, win=0.1,step=0.001)
      mat1 <- mat_ls[[2]]
      fig1 <- Heatmap(mat1, col=col_fun,
         cluster_rows=TRUE, cluster_columns=FALSE, row_km=4,
         show_row_names=FALSE, show_column_names=FALSE,
         column_title=paste("pseudotime(", treat1, ")", sep=""),
         heatmap_legend_param=list(title=""),
         use_raster=TRUE, raster_device="png")

     figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
                    "/Figure3_", oneMCl, "_lda", lda, "_trt", treat1, ".Win10Step0.1",
                    ".cluster.png", sep="")
     png(figfn, height=800, width=600, res=150)
     fig1 <- draw(fig1)
     r.list <- row_order(fig1)
     r.dend <- row_dend(fig1) 
     dev.off()
     ##
     clu_df <- lapply(names(r.list),function(i){
       out <- data.frame(genes=rownames(mat1)[r.list[[i]]],
                         Cluster=i, stringsAsFactors=FALSE)
       out                  
     })%>%do.call(rbind,.)
     rownames(clu_df) <- clu_df$genes
     #newCl <- c("1", "2", "3", "4")
     clu_df$New_cluster <- newCl[clu_df$Cluster] 
     opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
        "/Figure3_", oneMCl, "_lda", lda, "_trt", treat1, ".Win10Step0.1",
        ".cluster.csv", sep="")
     write.csv(clu_df, opfn)
###
### heatmap for condition that was contrated and keep row the same to that above

      treat0 <- contrast[1]
      mat0 <- mat_ls[[1]]
      ###
      mat0 <- mat0[clu_df$genes,]
      fig2 <- Heatmap(mat0, col=col_fun,
         cluster_rows=FALSE, cluster_columns=FALSE,
         show_row_names=FALSE, show_column_names=FALSE,
         row_order=clu_df$genes, row_title="gene",
         #row_split=factor(clu_df$Cluster, levels=1:4),
         #cluster_row_slices=FALSE,
         column_title=paste("pseudotime(", treat0, ")", sep=""),
         heatmap_legend_param=list(title=""),
         use_raster=FALSE, raster_device="png")

      figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/Figure3_",
              oneMCl, "_lda", lda, "_trt", treat0, ".Win10Step0.1", ".png", sep="")
         png(figfn, height=800, width=550, res=150)
         draw(fig2)
         dev.off()

###cluster
     ## fig3 <- Heatmap(mat0, col=col_fun,
     ##     cluster_rows=TRUE, cluster_columns=FALSE,
     ##     show_row_names=FALSE, show_column_names=FALSE,
     ##     row_title="gene",
     ##     column_title=paste("pseudotime(LDA_", lda, ")", sep=""),
     ##     heatmap_legend_param=list(title=""),
     ##     use_raster=TRUE, raster_device="png")

     ## figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
     ##                "/Figure3_", oneMCl, "_lda", lda, "_trt", treat0,
     ##                ".cluster.png", sep="")
     ## png(figfn, height=800, width=700, res=150)
     ## draw(fig3)
     ## dev.off() 
  } ##end for lda
} ## end for cell types    



#########
### 7 ###
#########
##
if(FALSE){
getData <- function(MCls, index=1, top=50){

   df2 <- map_dfr(MCls, function(oneMCl){
      fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/4.2_LDA.", oneMCl, ".loading.rds", sep="")
      w <- read_rds(fn)
      wi <- sort(abs(w[,index]), decreasing=T)
      gene0 <- gsub("\\.[0-9]*", "", names(wi))
      gene0 <- gene0[1:top]
      df0 <- bitr(gene0, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
      df0 <- df0%>%
             mutate(MCls=oneMCl)%>%
             filter(!grepl("LOC",SYMBOL))%>%distinct(ENSEMBL,.keep_all=T)
      df0$MCls <- oneMCl
      df0
   })   
   df2
}

###get data
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- getData(MCls, index=2, top=100)
opfn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/7.2.LDA2_top100.rds"
write_rds(df2, opfn)

### enrich results
load("./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData")
geneBG <- gsub("\\.[0-9].*", "", rownames(YtX))
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


fn <- "./9_RNA.dynamic2_output/Old/7.2.1_top100.rds"
df2 <- read_rds(fn)
cg <- compareCluster(ENTREZID~MCls, 
                     data=df2, 
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     universe=BgDf$ENTREZID, 
                     ont="ALL",
                     pvalueCutoff=0.1,
                     qvalueCutoff=0.5,
                     minGSSize=0,
                     maxGSSize=nrow(BgDf))
opfn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/7.2.2_enrichGO.rds"
write_rds(cg,opfn)

### plots
cg <- read_rds("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/7.2.2_enrichGO.rds")                 
p1 <- dotplot(cg, x=~MCls, showCategory=5)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))

figfn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/Figure8.2.1.png"
png(figfn,width=1000, height=1000, res=150)
print(p1)
dev.off() 

cg2 <- cg%>%
       filter(grepl("glucocorticoid|corticosteroid|lipopolysaccharide", Description))
p2 <- dotplot(cg2, x=~MCls, showCategory=5)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))

figfn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/Figure8.2.2.png"
png(figfn,width=1300, height=1000, res=150)
print(p2)
dev.off() 

### number of cells
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell") 
#prefix <- "./9_RNA.dynamic2_output/Old/LDA2Bin/"
#for (oneMCl in MCls){
#fn <- paste("./9_RNA.dynamic2_output/Old/4_LDA.", oneMCl, ".meta.rds", sep="")
#meta0 <- read_rds(fn)
#bti <- meta0 %>% transmute(bti=paste(BEST.GUESS, MCls, treats, BATCH, Bin2, sep="_")) %>% unlist %>% factor
#ncell <- data.frame(x=bti)%>%
#         group_by(x)%>%
#         summarise(ncell=n())
#opfn <- paste(prefix, "ncell.", oneMCl, ".RData", sep="")
#save(ncell,file=opfn)
#}
} ##

###
###
