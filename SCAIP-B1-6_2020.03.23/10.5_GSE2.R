###
library(tidyverse)
library(purrr)
library(qvalue)
##
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
theme_set(theme_grey())


rm(list=ls())

outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)



prefix <- "4_mu"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/"
 
### background gene list
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### geneCluster for enrichment results
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0, "Up", "Down"))



cg <- compareCluster(ENTREZID~contrast+MCls,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- paste(outdir, prefix, "_enrichGO.rds", sep="")
write_rds(cg, opfn1)
                                                 
### (4)
### KEGG enrichment analysis
ck <- compareCluster(ENTREZID~contrast+MCls,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     minGSSize=0,
                     maxGSSize=1000)
ck <- ck%>%mutate(Cluster1=Cluster)
opfn2 <- paste(outdir, prefix, "_enrichKEGG.rds", sep="")
write_rds(ck, opfn2)



####################
### scatterplots ###
####################
rm(list=ls())
cg <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/3_phiNew_enrichGO.rds")
cg <- cg%>%as.data.frame()%>%
   mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
          Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
          Diff.not=Diff.total-Diff.in,
          Bg.in=as.numeric(gsub("/.*","", BgRatio)),
          Bg.total=as.numeric(gsub(".*/","", BgRatio)),
          Bg.not=Bg.total-Bg.in,
          maxSize=as.numeric(gsub("/.*","", BgRatio)))%>%
   filter(maxSize>5&maxSize<500)
##
cg_mu <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/4_mu_enrichGO.rds")
cg_mu <- cg_mu%>%as.data.frame()%>%
   mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
          Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
          Diff.not=Diff.total-Diff.in,
          Bg.in=as.numeric(gsub("/.*","", BgRatio)),
          Bg.total=as.numeric(gsub(".*/","", BgRatio)),
          Bg.not=Bg.total-Bg.in,
          maxSize=as.numeric(gsub("/.*","", BgRatio)))%>%
   filter(maxSize>5&maxSize<500)
###
## shared <- c("type I interferon signaling pathway",
##    "response to lipopolysaccharide",
##    "cytokine-mediated signaling pathway",
##    "innate immune response")

##
cg2 <- cg%>%
   mutate(rn=paste(MCls, contrast, ID, sep="_"))%>%
   dplyr::select(rn, contrast, MCls, ID, Description, pvalue, p.adjust) 
###
mu2 <- cg_mu%>%
   mutate(rn=paste(MCls, contrast, ID, sep="_"))%>%
   dplyr::select(rn, pvalue, p.adjust)
df <- cg2%>%inner_join(mu2, by="rn")
df <- df%>%mutate(x=-log10(pvalue.x), y=-log10(pvalue.y))

###
## df <- df%>%mutate(gr=ifelse(Description%in%shared, 1, 2))

p <- ggplot(df, aes(x=x, y=y))+
   rasterise(geom_point(colour="grey30", size=0.5),dpi=300)+ 
   ## geom_point(aes(colour=factor(gr), size=factor(gr)))+
   ## scale_colour_manual(values=c("1"="green", "2"="grey30"))+
   ## scale_size_manual(values=c("1"=1.5, "2"=0.5))+ 
   geom_abline(colour="red")+ 
   facet_grid(MCls~contrast,scales="free",
      labeller=labeller(contrast=c("LPS"="LPS", "LPS-DEX"="LPS+DEX","PHA"="PHA", "PHA-DEX"="PHA+DEX")))+ 
   xlab(bquote(~log[10]~"(pvalue)"~" enrich for DVG"))+
   ylab(bquote(~log[10]~"(pvalue)"~" enrich for DEG"))+
   theme_bw()+
   theme(## legend.position="none",
         axis.title=element_text(size=10),
         axis.text=element_text(size=8),
         strip.text=element_text(size=12))

## figfn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/Figure1.1_GO.scatter.png"
## png(figfn, width=750, height=750, res=120)
## print(p)
## dev.off()

figfn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler2/Figure1.1_GO.scatter.pdf"
pdf(figfn, width=8, height=8)
print(p)
dev.off()
