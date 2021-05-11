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
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F)


############################################################
### 2021-05-03, psedotime definition based on 6,571 DEGs ###
###       diagonal linear Discriminant analysis          ###
###       Last modified by Julong wei, 2021-05-03        ###
############################################################

rm(list=ls())

########################
### 1. calculate LDA ### 
########################

res <- read_rds("./6_DEG.CelltypeNew_output/Filter2/3_Batch1456.meta.rds")
i <- 1
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
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
                       "LPS_DEX"=c("LPS-EtOH", "LPS-DEX"),
                       "PHA"=c("CTRL", "PHA-EtOH"),
                       "PHA_DEX"=c("PHA-EtOH", "PHA-DEX"))
   contrast_nn <- names(contrast_ls)
   metaNew <- sapply(contrast_nn, function(nn){
      oneX <- contrast_ls[[nn]] 
      ### extract lfc as weights
      nn2 <- gsub("_","-",nn) ### change LPS_DEX into LPS+DEX 
      res0 <- res%>%filter(MCls==oneMCl, contrast==nn2)
      lfc <- abs(res0$beta)
      lfc[is.na(lfc)] <- 0
      names(lfc) <- as.character(res0$gene)
      lfc0 <- lfc[rn]
      
      ### matrix 1 
      ii <- meta%>%filter(treats==oneX[1])%>%dplyr::pull(NEW_BARCODE)
      X1 <-X[,ii]
      mu1 <- apply(X1, 1, mean)

      ### matrix 2
      ii <- meta%>%filter(treats==oneX[2])%>%dplyr::pull(NEW_BARCODE)
      X2 <-X[,ii]
      mu2 <- apply(X2, 1, mean)

      Diff <- as.matrix((mu2-mu1)*(1/var)) ##
      Diff <- as.matrix((mu2-mu1)*(1/var)*lfc0) ##shrinkage
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
   opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/2_MCls.",
                  oneMCl, ".SDLDA.rds", sep="")
   write_rds(metaNew, opfn)
      
   s2 <- Sys.time()
   d12 <- difftime(s2, s1, units="mins")
   cat(oneMCl, ":", d12, "\n")
}, mc.cores=1)



#########################
### summary LDA plots ###
#########################
### scatter plots plus density plots


col1 <- c("CTRL"="#828282",
          "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
          "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
df2 <- map_dfr(MCls, function(oneMCl){
   ##
   cat(oneMCl, "\n")
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/2_MCls.",
               oneMCl, ".SDLDA.rds", sep="")
   dd <- read_rds(fn)%>%mutate(MCls=oneMCl)
   dd
})

### (1)., LPS vs PHA
## together
fig_ls <- lapply(MCls, function(oneMCl){
   df0 <- df2%>%filter(MCls==oneMCl)
   fig0 <- ggplot(df0, aes(x=z_LPS, y=z_PHA))+
      geom_point(aes(colour=factor(treat2)), size=0.1)+
      scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
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
   ggplot(df2%>%filter(MCls=="Bcell"), aes(z_LPS, z_PHA))+
   geom_point(aes(colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_rect(colour=NA, fill=NA),
         legend.text=element_text(size=8),
         legend.key=element_rect(fill=NA),
         legend.key.size=grid::unit(1,"lines"))
   )

##
png("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/Figure2.X1.png", width=900, height=700, res=130)

fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


### (2)., LPS-DEX vs PHA-DEX
## together
fig_ls <- lapply(MCls, function(oneMCl){
   df0 <- df2%>%filter(MCls==oneMCl)
   fig0 <- ggplot(df0, aes(x=z_LPS_DEX, y=z_PHA_DEX))+
      geom_point(aes(colour=factor(treat2)), size=0.1)+
      scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
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
   ggplot(df2%>%filter(MCls=="Bcell"), aes(z_LPS_DEX, z_PHA_DEX))+
   geom_point(aes(colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_rect(colour=NA, fill=NA),
         legend.text=element_text(size=8),
         legend.key=element_rect(fill=NA),
         legend.key.size=grid::unit(1,"lines"))
   )

##
png("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/Figure2.X2.png", width=900, height=700, res=130)

fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()



### (3)., LPS vs LPS-DEX
## together
fig_ls <- lapply(MCls, function(oneMCl){
   df0 <- df2%>%filter(MCls==oneMCl)
   fig0 <- ggplot(df0, aes(x=z_LPS, y=z_LPS_DEX))+
      geom_point(aes(colour=factor(treat2)), size=0.1)+
      scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
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
   ggplot(df2%>%filter(MCls=="Bcell"), aes(z_LPS, z_LPS_DEX))+
   geom_point(aes(colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_rect(colour=NA, fill=NA),
         legend.text=element_text(size=8),
         legend.key=element_rect(fill=NA),
         legend.key.size=grid::unit(1,"lines"))
   )

##
png("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/Figure2.X3.png", width=900, height=700, res=130)

fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()


### (4)., PHA vs PHA-DEX
## together
fig_ls <- lapply(MCls, function(oneMCl){
   df0 <- df2%>%filter(MCls==oneMCl)
   fig0 <- ggplot(df0, aes(x=z_PHA, y=z_PHA_DEX))+
      geom_point(aes(colour=factor(treat2)), size=0.1)+
      scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=2)))+
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
   ggplot(df2%>%filter(MCls=="Bcell"), aes(z_PHA, z_PHA_DEX))+
   geom_point(aes(colour=factor(treat2)), size=0.1)+
   scale_colour_manual(values=col1, guide=guide_legend(override.aes=list(size=1)),labels=lab2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_rect(colour=NA, fill=NA),
         legend.text=element_text(size=8),
         legend.key=element_rect(fill=NA),
         legend.key.size=grid::unit(1,"lines"))
   )

##
png("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/Figure2.X4.png", width=900, height=700, res=130)

fig1 <- plot_grid(fig_ls[[1]], fig_ls[[2]],
                  fig_ls[[3]], fig_ls[[4]],
                  nrow=2, ncol=2, align="hv",axis="tb")

print(plot_grid(fig1, legend2, rel_widths=c(4,1)))
dev.off()



################################################
### 2. defined Three bins according to 4 LDAs ###
################################################

rm(list=ls())

outdir <- "./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/DLDA_Bin/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F)

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
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/1_MCls.", oneMCl, ".DLDA.rds", sep="")
   meta <- read_rds(fn)
   treats <- sort(unique(meta$treats)) 

   meta <- map_dfr(treats, function(ii){
      meta0 <- meta%>%filter(treats==ii)
      meta0$Bin_LPS <- Binfun(meta0$z_LPS)
      meta0$Bin_LPS_DEX <- Binfun(meta0$z_LPS_DEX)
      meta0$Bin_PHA <- Binfun(meta0$z_PHA)
      meta0$Bin_PHA_DEX <- Binfun(meta0$z_PHA_DEX)
      meta0
   })
###
   opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/DLDA_Bin/1_MCls.", oneMCl, ".Bin.rds", sep="")
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
Bins <- c("LPS", "LPS.DEX", "PHA", "PHA.DEX")
###
for (oneMCl in MCls){

fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/DLDA_Bin/1_MCls.", oneMCl, ".Bin.rds", sep="")
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
   prefix <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/DiagLDA/DLDA_Bin/2_", oneMCl, "_Bin", i, ".", Bins[i], sep="")
   opfn1 <- paste(prefix,"_YtX.sum.rds", sep="")
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
