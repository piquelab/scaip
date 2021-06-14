#
###
library(tidyverse)
library(purrr)
library(furrr) 
library("BiocParallel")
library(qqman)
library(qvalue)
##
library(DESeq2)
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
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./7_GSE.ClusterProfiler_output/Filter2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#############################################################
### Enrichment analysis for differentially expressed gene ###
###     12-29-2020, last modified by Julong wei,           ###
#############################################################


##################################################################
### Example code for enrichment analysis using ClusterProfiler ###
##################################################################
#load("./6_DEG_CelltypeNew_output/YtX.comb.RData")
#geneBG <- rownames(YtX)
#geneBG.df <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

###example
#load("./6_DEG_CelltypeNew_output/Sigs.gene.DEG.RData")
#gene <- sigs[["LPS"]][["Monocyte"]]
#gene.df <- bitr(gene,fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"),OrgDb=org.Hs.eg.db)
  
#ego2 <- enrichGO(gene=gene.df$ENTREZID,
#                 universe=geneBG.df$ENTREZID,
#                 #keyType="ENSEMBL",
#                 OrgDb=org.Hs.eg.db,
#                 ont="ALL",
#                 minGSSize=1,
#                 maxGSSize=nrow(gene.df),
#                 pvalueCutoff=1, qvalueCutoff=1, readable=T)
#ego3 <- gseGO(gene=gene.df$ENTREZID,
#              universe=geneBG.df$ENTREZID,
#              OrgDb=org.Hs.eg.db,
#              ont="ALL",
#              nPerm=1000,
#              minGSSize=1,
#              maxGSSize=nrow(gene.df),
#              pvalueCutoff=1, qvalueCutoff=1, verbose=FALSE)


#####################################################
### 0, parse gtt3 and obtain protein coding genes ###
#####################################################
#if(F){
#rm(list=ls())
#cat("0.", "Extract protein coding genes", "\n")
##/wsu/home/groups/piquelab/SCAIP/covariates/gencode.v34lift37.basic.annotation.gff3.gz
#fn <- "/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz"
#con <- gzfile(fn,"rt")
#gff3 <- read.table(con,header=F)
#ID <- gsub(".*gene_id=", "", gff3$V9)
#ID <- gsub(";.*", "", ID)
#gene_type <- gsub(".*gene_type=", "", gff3$V9)
#gene_type <- gsub(";.*","",gene_type)
#
#anno <- gff3[,-9]%>%
#        mutate(ID=ID, gene_type=gene_type)%>%
#        filter( grepl("protein_coding",gene_type))%>%
#        distinct(ID, .keep_all=T)
#
#### 19,957 protein coding gene                 
#pList <- gsub("\\..*", "", anno$ID) 
#}

#######################################
### 1, gene set enrichment analysis ###
#######################################
#if (FALSE){
#### background gene list
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#geneBG <- gsub("\\..*", "", rownames(YtX))
#geneBG0 <- geneBG[geneBG%in%pList]  ##17,690
#geneBG.df <- bitr(geneBG0, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
#
###DEG gene list, Among 6413 DEG genes, 5930 protein gene
#load("./6_DEG_CelltypeNew_output/Sigs.gene.DEG.RData")
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
#geneList <- list()
#for(icont in contrast){
#   for (oneMCl in MCls){
#      gene <- sigs[[icont]][[oneMCl]]
#      gene0 <- gene[gene%in%pList]
#      gene.df <- bitr(gene0, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
#      comb <- paste(icont, oneMCl, sep="_")      
#      geneList[[comb]] <- gene.df$ENTREZID
#   }
#}
#   
#
#cg <- compareCluster(geneCluster=geneList, 
#                     universe=geneBG.df$ENTREZID,
#                     fun="enrichGO", 
#                     OrgDb="org.Hs.eg.db", 
#                     ont="ALL")
#                     
#fig0 <- dotplot(cg%>%filter(ONTOLOGY=="BP"),showCategory=5)+
#        theme(axis.text.x=element_text(angle=90, hjust=1,size=10),
#              axis.text.y=element_text(size=10))
#
#figfn <- "./7_GSE.ClusterProfiler_output/Figure1.1.BP.png"
#png(figfn, width=1200, height=1000, res=120)
#fig0
#dev.off()
#}


######################################################
### 1, Enrichment for up and down genes separately ###
######################################################

### background gene list
load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
rownames(YtX_sel) <- gsub("\\.[0-9].*", "", rownames(YtX_sel))

grch38_unq <- grch38%>%distinct(ensgene, .keep_all=T)
anno <- data.frame(ensgene=rownames(YtX_sel))%>%mutate(rnz=rowSums(YtX_sel))%>%left_join(grch38_unq, by="ensgene")
geneBG <- anno%>%dplyr::pull(ensgene)##15,770
#geneBG0 <- geneBG[geneBG%in%pList]  ##17,690
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)


fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0,"Up", "Down"))


### GO enrichment 
# "GO analysis for direction (up and down) directly", "\n")
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

#cg0 <- as.data.frame(cg)
cg <- cg%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(cg,"./7_GSE.ClusterProfiler_output/Filter2/old/1_enrichGO.rds")

###
### combine direction of up and down
## "GO analysis for direction together" 
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

write_rds(cg,"./7_GSE.ClusterProfiler_output/Filter2/1.1_enrichGO.combDirection.rds")

                                                 
###(3), KEGG enrichment analysis
cat("(3)", "KEGG analysis", "\n")
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=1,
                     qvalueCutoff=1,
                     minGSSize=0,
                     maxGSSize=1000)

#ck0 <- as.data.frame(ck)
ck <- ck%>%mutate(Cluster1=gsub("\\.((Down)|(Up))", "", Cluster))
write_rds(ck,"./7_GSE.ClusterProfiler_output/Filter2/old/2_enrichKEGG.rds")
              

######################################################
### 2. Show general figures of enrichment results  ###
######################################################

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
cluster2 <- setNames(cl2, cl)
lab2 <- setNames(gsub("-","+",cl),cl2)
        
                      
###(1), showing all the GO
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)]  
cg2 <- cg%>%mutate(ClusterNew=x,
                   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
   filter(maxGSSize<500, p.adjust<0.1)                  
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))

### pdf
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.GO.pdf"
pdf(figfn, width=15, height=10)
print(fig1)
dev.off()

### png
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure1.GO.png"
png(figfn, width=3500, height=2000, res=180)
print(fig1)
dev.off()


### (3). showing KEGG
ck <- read_rds("7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x,
                   maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
  filter(maxGSSize<500,p.adjust<0.1)

fig3 <- enrichplot::dotplot(ck2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=10))
        
### pdf
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure3.KEGG.pdf"
pdf(figfn, width=15, height=10)
print(fig3)
dev.off()

### png
figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure3.KEGG.png"
png(figfn, width=3000, height=1500, res=150)
print(fig3)
dev.off()

##################################
### 3. show specific GO terms  ###
##################################


### background gene list
#load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
#rownames(YtX_sel) <- gsub("\\.[0-9].*", "", rownames(YtX_sel))
#geneBG <- rownames(YtX_sel)
#BG <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
#
#cg <- enrichGO(gene=BG$ENTREZID,
#                 universe=BG$ENTREZID,
#                 OrgDb=org.Hs.eg.db,
#                 ont="ALL",
#                 minGSSize=0,
#                 maxGSSize=nrow(BG),
#                 pvalueCutoff=1, qvalueCutoff=1)
#write_rds(cg, file="./7_GSE.ClusterProfiler_output/Filter2/4_GO.gene.rds")


#gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
#ens <- gene2%>%dplyr::pull(ENSEMBL)

ExampleGOplot <- function(cg){
   
   x <- str_split(cg$GeneRatio, "/", simplify=T)
   GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
   Drt2 <- c("Up"=1, "Down"=2) 
   cg <- cg%>%mutate(Direction2=Drt2[direction], 
              contrast2=paste(contrast, Direction2, sep="."),
              contrast2=gsub("-", "+", contrast2))
   #             
   cg$size <- rep(1,nrow(cg))
   cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
   cg$size[GeneRatio>=0.15] <- 3
   #
   ## cg$p2 <- cg$p.adjust
   ## cg$p2[cg$p2>0.1] <- NA
   
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls))+
           geom_point(aes(size=factor(size), colour=p2))+
           scale_x_discrete(labels=c("LPS.1"="LPS.Up", "LPS.2"="LPS.Down",
                            "LPS+DEX.1"="LPS+DEX.Up", "LPS+DEX.2"="LPS+DEX.Down",
                            "PHA.1"="PHA.Up", "PHA.2"="PHA.Down",
                            "PHA+DEX.1"="PHA+DEX.Up", "PHA+DEX.2"="PHA+DEX.Down"))+
           scale_colour_gradient(name="p.adjust",                           
                                 low="blue", high="red", trans="reverse")+    #"#ffa500"
           scale_size_manual(name="GeneRatio", 
                             values=c("1"=1, "2"=2, "3"=3),
                             labels=c("1"="<0.05", "2"="<0.15", "3"=">=0.15"))+
           theme_bw()+
           theme(axis.title=element_blank(),
                 axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
                 axis.text.y=element_text(size=8),
                 legend.background=element_blank(),
                 legend.title=element_text(size=8),
                 legend.text=element_text(size=6),
                 legend.key.size=grid::unit(0.8,"lines"))
   fig0
}

### type I interferon
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.1] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("Type I interferon signaling pathway")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.1_interferon.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 

### response to glucocorticoid ###
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("cellular response to glucocorticoid stimulus", Description))
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("cellular response to glucocorticoid stimulus")+
        theme(plot.title=element_text(hjust=0.5, size=8))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.2_glucocorticoid.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


### inflammatory response
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="inflammatory response")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("inflammatory response")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.3_inflammatory.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off()


### cell activation involved in immune response 
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cell activation involved in immune response")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("Cell activation involved in immune response")+
        theme(plot.title=element_text(hjust=0.5, size=8))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.4_immuneresponse.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off()


### innate immune response
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="innate immune response")
##
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("innate immune response")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.5_immuneresponse.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


### response to bacterium 

cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="response to bacterium")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("response to bacterium")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.6_bacterium.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


### cytokine-mediated signaling pathway
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine-mediated signaling pathway")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("cytokine-mediated signaling pathway")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.7_cytokine-mediated.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


### cytokine receptor binding
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine receptor binding")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("cytokine receptor binding")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.8_cytokine-receptor.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


###
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(Description=="cytokine activity")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.05] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("cytokine actvity")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.9_cytokine-activity.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


###
### response to lipopolysaccharide
cg2 <- cg0%>%filter(Description=="response to lipopolysaccharide")
##
cg2$p2 <- cg2$p.adjust
cg2$p2[cg2$p2>0.1] <- NA

fig1 <- ExampleGOplot(cg2)+
        ggtitle("response to lipopolysaccharide")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure4.10_lipopoly.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 


############
### KEGG ###
############

ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds") 
ck0 <- as.data.frame(ck)

### Influenza A
ck2 <- ck0%>%filter(Description=="Influenza A")
ck2$p2 <- ck2$p.adjust
ck2$p2[ck2$p2>0.1] <- NA

fig1 <- ExampleGOplot(ck2)+
        ggtitle("Influenza A")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.1_Influenza_A.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off() 

### cytokine-cytokine 
ck2 <- ck0%>%filter(Description=="Cytokine-cytokine receptor interaction")
ck2$p2 <- ck2$p.adjust
ck2$p2[ck2$p2>0.1] <- NA

fig1 <- ExampleGOplot(ck2)+
        ggtitle("Cytokine-cytokine receptor interaction")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.2_Cytokine.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off()

### COVID-19
ck2 <- ck0%>%filter(Description=="Coronavirus disease - COVID-19")
ck2$p2 <- ck2$p.adjust
ck2$p2[ck2$p2>0.1] <- NA

fig1 <- ExampleGOplot(ck2)+
        ggtitle("Coronavirus disease-COVID-19")+
        theme(plot.title=element_text(hjust=0.5, size=10))

figfn <- "./7_GSE.ClusterProfiler_output/Filter2/Figure5.3_COVID-19.png"
png(figfn, width=500, height=400, res=120)
print(fig1)
dev.off()

###
#avePathway <- function(X){
#   bti <- colnames(X)
#   cvt <- str_split(bti, "_", simplify=T)
#   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])%>%
#          mutate(comb=paste(MCls, Batch, sep="_"))
#   comb <- unique(cvt$comb)
#   for (ii in comb){
#      bti0 <- cvt%>%filter(comb==ii)%>%dplyr::pull(rn)
#      x <- X[,bti0]
#      x.mean <- apply(x, 1, mean, na.rm=T)
#      x.scale <- sweep(x, 1, x.mean, "-")
#      X[,bti0] <- x.scale
#   }
#   ###
#   MCls <- unique(cvt$MCls)
#   X2 <- map_dfc(MCls, function(ii){
#      ###LPS and PHA
#      bti0 <- cvt%>%filter(MCls==ii,treats=="CTRL")%>%dplyr::pull(rn)
#      x0 <- apply(x[,bti0], mean, na.rm=T)
#      bti0 <- cvt%>%filter(MCls==ii,treats=="LPS-EtOH")%>%dplyr::pull(rn)
#      xLPS <- sweep(x[,bti0], 1, x0, "-")
#      bti0 <- cvt%>%filter(MCls==ii,treats=="PHA-EtOH")%>%dplyr::pull(rn)
#      xPHA <- sweep(x[,bti0], 1, x0, "-")
#      ###LPS-DEX
#      bti0 <- cvt%>%filter(MCls==ii, treats=="LPS-EtOH")%>%dplyr::pull(rn)
#      x0 <- apply(x[,bti0], 1, mean, na.rm=T)
#      bti0 <- cvt%>%filter(MCls==ii,treats=="LPS-DEX")%>%dplyr::pull(rn)
#      xLPS-DEX <- sweep(x[,bti0], 1, x0, "-")
#      ###PHA-DEX
#      bti0 <- cvt%>%filter(MCls==ii, treats=="PHA-EtOH")%>%dplyr::pull(rn)
#      x0 <- apply(x[,bti0], 1, mean, na.rm=T)
#      bti0 <- cvt%>%filter(MCls==ii,treats=="PHA-DEX")%>%dplyr::pull(rn)
#      xPHA-DEX <- sweep(x[,bti0], 1, x0, "-")
#      ###
#      x <- cbind(xLPS,xLPS-DEX, xPHA, xPHA-DEX)                  
#    })  
#       
#   
#   pathway <- apply(X, 2, mean, na.rm=T)
#}
#
####bulk, NB.mu, NB.phi
#getPathwayData <- function(gene, datatype="bulk"){
#
#### bulk data
#if (datatype=="bulk"){
#   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
#   load(fn)
#   rn <- gsub("\\..*", "", rownames(YtX_sel))
#   rownames(YtX_sel) <- rn
#   X <- YtX_sel[rn%in%gene,]+1
#
#   bti <- colnames(YtX_sel)
#   cvt <- str_split(bti, "_", simplify=T)
#   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])
#
#   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
#   load(fn)
#   counts <- colSums(YtX)
#   counts <- counts[colnames(YtX_sel)]
#
#   X <- sweep(X, 2, counts,"/")
#   X <- X*1e+06 
#   X <- log2(X)
#   cvt$y <- avePathway(X)
#} ###
#
#### NB.mu
#if(datatype=="NB.mu"){
####
#   load("./10_RNA.Variance_output/tmp10/1.2_Sel.Bx.RData")
#   rn <- gsub("\\..*", "", rownames(Bx))
#   rownames(Bx) <- rn
#
#   X <- Bx[rn%in%gene,]
#
#   bti <- colnames(Bx)
#   cvt <- str_split(bti, "_", simplify=T)
#   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])
#   X <- log2(X)
#   cvt$y <- avePathway(X) 
#}
#
#### NB.phi
#if(datatype=="NB.phi"){
####
#   load("./10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData")
#   rn <- gsub("\\..*", "", rownames(PhxNew2))
#   rownames(PhxNew2) <- rn
#   X <- PhxNew2[rn%in%gene,]
#
#   bti <- colnames(PhxNew2)
#   cvt <- str_split(bti, "_", simplify=T)
#   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=cvt[,2], sampleID=cvt[,3], Batch=cvt[,4])
#   X <- log2(X)
#   cvt$y <- avePathway(X) 
#}
#cvt
#}


