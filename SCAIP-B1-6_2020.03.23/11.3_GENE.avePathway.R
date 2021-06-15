###
###
library("rhdf5")
library("corpcor")
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(purrr)
##
library(Seurat)
library(SeuratDisk)
library(annotables)
library(biobroom)
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())


###############################################
#### Average expression in given pathway    ###
###############################################
### last modified 6-14-2021, Julong wei
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#geneBG <- gsub("\\.[0-9].*", "", rownames(YtX))
#Bg <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

##
avePathway <- function(X){

### filtering more missing value   
   ii <- apply(!is.na(X), 1, sum)
   X <- X[ii>20,]
   bti <- colnames(X)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1],treats=cvt[,2],sampleID=cvt[,3],Batch=cvt[,4])%>%
          mutate(comb=paste(MCls, Batch, sep="_"))
   comb <- unique(cvt$comb)
   for (ii in comb){
      bti0 <- cvt%>%filter(comb==ii)%>%dplyr::pull(rn)
      x <- X[,bti0]
      x.mean <- apply(x, 1, mean, na.rm=T)
      x.scale <- sweep(x, 1, x.mean, "-")
      X[,bti0] <- x.scale
   }
   pathway <- apply(X, 2, mean, na.rm=T)
}

###bulk, NB.mu, NB.phi
getPathwayData <- function(gene, datatype="bulk"){

### bulk data
if (datatype=="bulk"){
   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
   load(fn)
   rn <- gsub("\\..*", "", rownames(YtX_sel))
   rownames(YtX_sel) <- rn
   X <- YtX_sel[rn%in%gene,]+1

   bti <- colnames(YtX_sel)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2],sampleID=cvt[,3],Batch=cvt[,4])

   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
   load(fn)
   counts <- colSums(YtX)
   counts <- counts[colnames(YtX_sel)]

   X <- sweep(X, 2, counts,"/")
   X <- X*1e+06 
   X <- log2(X)
   cvt$y <- avePathway(X)
} ###

### NB.mu
if(datatype=="NB.mu"){
###
   load("./10_RNA.Variance_output/tmp10/1.2_Sel.Bx.RData")
   rn <- gsub("\\..*", "", rownames(Bx))
   rownames(Bx) <- rn

   X <- Bx[rn%in%gene,]

   bti <- colnames(Bx)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])
   X <- log2(X)
   cvt$y <- avePathway(X) 
}

### NB.phi
if(datatype=="NB.phi"){
###
   load("./10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData")
   rn <- gsub("\\..*", "", rownames(PhxNew2))
   rownames(PhxNew2) <- rn
   X <- PhxNew2[rn%in%gene,]

   bti <- colnames(PhxNew2)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])
   X <- log2(X)
   cvt$y <- avePathway(X) 
}
cvt
}


### 1, type I interferon signaling pathway
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Type I interferon signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure1.1_typeIinterferon.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Type I interferon signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure1.2_typeIinterferon.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Type I interferon  signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure1.3_typeIinterferon.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        

###(2), cytokine receptor binding
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine receptor binding")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine receptor binding score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.1_cytokineReceptorBinding.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine receptor binding score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.2_cytokineReceptorBinding.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine receptor binding score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.3_cytokineReceptorBinding.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()        


###
### (3), cytokine-mediated signaling pathway
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine-mediated signaling pathway")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine-mediated signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.1_cytokine-mediated.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine-mediated signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.2_cytokine-mediated.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine-mediated signaling pathway score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.3_cytokine-mediated.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()


### (4), cytokine activity
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine activity")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
### bulk           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine activity score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure4.1_cytokineActivity.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()


### NB.mu
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine activity score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure4.2_cytokineActivity.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

### NB.phi
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("cytokine activity score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure4.3_cytokineActivity.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        

### (5), response to bacterium
outdir <- "./11_GENE.Example_output/6_avePathway/"
### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="response to bacterium")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

### (1)           
### bulk           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.1_responsebacterium.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

### (2)
### NB.mu
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.2_responsebacterium.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
### (3), NB.phi
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.3_responsebacterium.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()        


### (6).
### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="innate immune response")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("innate immune response score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.1_innateimmune.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("innate immune response score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.2_innateimmune.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("innate immune response score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.3_innateimmune.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        

### 7, response to lipopolysaccharide
### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="response to lipopolysaccharide")
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)

###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to lipopolysaccharide score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure7.1_lipopoly.bulk.png" 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to lipopolysaccharide score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure7.2_lipopoly.NB.mu.png"
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to lipopolysaccharide score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure7.3_lipopoly.NB.phi.png"
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()



###
###8, Influenza A
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds") 
ck0 <- as.data.frame(ck)
ck2 <- ck0%>%filter(Description=="Influenza A")
geneList <- unique(unlist(str_split(ck2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)


###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Influenza A score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure8.1_Influenza_A.bulk.png" 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Influenza A score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure8.2_Influenza_A.NB.mu.png"
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Influenza A score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure8.3_Influenza_A.NB.phi.png"
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()


###
###9, "Cytokine-cytokine receptor interaction"
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds") 
ck0 <- as.data.frame(ck)
ck2 <- ck0%>%filter(Description=="Cytokine-cytokine receptor interaction")
geneList <- unique(unlist(str_split(ck2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)


###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Cytokine-cytokine receptor interaction score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure9.1_Cytokine.bulk.png" 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Cytokine-cytokine receptor interaction score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure9.2_Cytokine.NB.mu.png"
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Cytokine-cytokine receptor interaction score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure9.3_Cytokine.NB.phi.png"
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()


###
###10, "Cytokine-cytokine receptor interaction"
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds") 
ck0 <- as.data.frame(ck)
ck2 <- ck0%>%filter(Description=="Coronavirus disease - COVID-19")
geneList <- unique(unlist(str_split(ck2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)


###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
cvt <- getPathwayData(ens, datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Coronavirus disease score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure10.1_COVID-19.bulk.png" 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Coronavirus disease score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure10.2_COVID-19.NB.mu.png"
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens, datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("Coronavirus disease score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- "./11_GENE.Example_output/6_avePathway/Figure10.3_COVID-19.NB.phi.png"
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
