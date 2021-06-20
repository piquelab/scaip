###
#library(Matrix)
#library(MASS)
library(tidyverse)
library(parallel)
library(data.table)
library("BiocParallel")
library(qqman)
library(qvalue)
##
library(annotables)
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
theme_set(theme_grey())


#####################
### enrich eGenes ###
#####################

rm(list=ls())

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
dataset <- data.frame(MCls=rep(MCls, each=4), LDA=rep(LDA, times=4))

contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                     "PHA"=c("CTRL", "PHA-EtOH"),
                     "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))


enrich.eGene <- function(df, eGene, n_tested=11980){
  n_eGene <- length(eGene)
  ncl <- sort(unique(df$Cluster))
  res <- map_dfr(ncl, function(i){
     g <- df%>%filter(Cluster==i)
     g1 <- g%>%filter(genes%in%eGene)    
     d <- data.frame(gene.in.interest=c(nrow(g1),nrow(g)-nrow(g1)),
       gene.not.interest=c(n_eGene, n_tested-n_eGene))
     rownames(d) <- c("In_category", "not_in_category")
     
    test <- fisher.test(d, alternative="greater")
    res <- data.frame(Cluster=i, n_gene=nrow(g), n_eGene=nrow(g1),
                      pvalue=test$p.value, OR=test$estimate) 
  })
  rownames(res)<- NULL
  res
}

###
eGene <- read_rds("/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/DiagLDA2/eGene_list.rds")
###
option <- "DiagLDA2"
for (i in c(13, 15, 16)){

  i <- 16  
  oneMCl <- dataset[i,1]
  lda <- dataset[i,2]
  contrast <- contrast_ls[[lda]]
  treat1 <- contrast[2]

  fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
    "/Cluster_genes/",     oneMCl, "_lda", lda, "_trt", treat1,
    ".cluster.csv", sep="")
  clu_df <- read.csv(fn)
  
  res_eGene_LDA <- enrich.eGene(clu_df, eGene[["eGene_LDA"]])

prop.test(x=c(13, 36), n=c(336, 413), alternative="less")
    

##############################
### enrichment analysis ###
##############################

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
LDA <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
dataset <- data.frame(MCls=rep(MCls, each=4), LDA=rep(LDA, times=4))
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                    "PHA"=c("CTRL", "PHA-EtOH"),
                    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
    
### background genes    
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"),OrgDb=org.Hs.eg.db)

###
option <- "DiagLDA2"
    
### Cluster genes
dfcomb <- map_dfr(13:16, function(i){
   oneMCl <- dataset[i,1]
   lda <- dataset[i,2]
   contrast <- contrast_ls[[lda]]
   treat1 <- contrast[2]
   fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
      "/Cluster_genes/", oneMCl, "_lda", lda, "_trt", treat1, ".cluster.csv", sep="")
   df <- read.csv(fn)%>%mutate(MCls=oneMCl, LDA=lda, treats=treat1)
})
    
### data directly used for enrichment analysis  
df0 <- bitr(dfcomb$genes, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- dfcomb%>%inner_join(df0,by=c("genes"="ENSEMBL"))%>%
  mutate(Cluster2=paste(treats, Cluster, sep="_"))

    
#####################
### GO enrichment ###
#####################

option <- "DiagLDA2"    
oneMCl <- "Tcell"
    
cg <- compareCluster(ENTREZID~Cluster2,
  data=geneCluster, universe=BgDf$ENTREZID,
  fun="enrichGO", OrgDb="org.Hs.eg.db",
  pvalueCutoff=1, qvalueCutoff=1, ont="ALL",
  minGSSize=0, maxGSSize=1000)

opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl,  "_enrichGO.rds", sep="")
write_rds(cg, opfn)    

###    
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl,  "_enrichGO.rds", sep="")    
cg <- read_rds(fn)    
cg2 <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
  ClusterNew=gsub("-EtOH", "", Cluster2))%>%
  filter(maxGSSize<500, p.adjust<0.1)

### png    
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=10))    
figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/",  oneMCl, "_Figure1.1_enrichGO.png", sep="")
png(figfn, width=1000, height=800, res=100)
print(fig1)
dev.off()    
    
### pdf    
fig2 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=10))
    
figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl,  "_Figure1.2_enrichGO.pdf", sep="")
pdf(figfn, width=12, height=10)
print(fig2)
dev.off()


    
################################
### KEGG enrichment analysis ###
################################
    
ck <- compareCluster(ENTREZID~Cluster2,
  data=geneCluster, universe=BgDf$ENTREZID,
  fun="enrichKEGG",
  pvalueCutoff=1, qvalueCutoff=1,
  minGSSize=0, maxGSSize=1000)

opfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl, "_enrichKEGG.rds", sep="")
write_rds(ck, opfn)    

fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl, "_enrichKEGG.rds", sep="")
ck <- read_rds(fn)    
ck2 <- ck%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)),
  ClusterNew=gsub("-EtOH", "", Cluster2))%>%
  filter(maxGSSize<500, p.adjust<0.1)    
    
### png
fig1 <- enrichplot::dotplot(ck2, x=~ClusterNew, showCategory=5)+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=10))    
figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl, "_Figure1.1_enrichKEGG.png", sep="")
png(figfn, width=1000, height=800, res=100)
print(fig1)
dev.off()


### pdf
fig2 <- enrichplot::dotplot(ck2, x=~ClusterNew, showCategory=10)+
  theme(axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=10))    
figfn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/", option,
  "/Cluster_genes/", oneMCl, "_Figure1.2_enrichKEGG.pdf", sep="")
pdf(figfn, width=12, height=10)
print(fig2)
dev.off()    


