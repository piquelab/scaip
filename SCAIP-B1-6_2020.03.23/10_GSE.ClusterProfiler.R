rm(list=ls())
source("./Bin/LibraryPackage.R")
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

######################################################
### 1, Enrichment for up and down genes separately ###
######################################################
if (FALSE){

cat("2.", "Enrichment analysis", "\n")
prefix <- "3_phiNew"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"

### (1) 
### background gene list
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)
geneBG <- unique(res$gene)
BgDf <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### (2)
### geneCluster for enrichment results
fn <- paste("./10_RNA.Variance_output/tmp9/", prefix, ".meta", sep="")
res <- read.table(fn, header=T)%>%drop_na(beta,qval)%>%filter(abs(beta)>0.5, qval<0.1)
df0 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)
geneCluster <- res%>%inner_join(df0,by=c("gene"="ENSEMBL"))%>%mutate(direction=ifelse(beta>0, "Up", "Down"))

### (3) 
### GO enrichment 
cat("2.1", "GO analysis", "\n")
cg <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=0.2,
                     qvalueCutoff=0.5, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- paste(outdir, prefix, ".enrichGO.2.rds", sep="")
write_rds(cg, opfn1)
                                                 
### (4)
### KEGG enrichment analysis
cat("2.2", "KEGG analysis", "\n")
ck <- compareCluster(ENTREZID~contrast+MCls+direction,
                     data=geneCluster,  
                     universe=BgDf$ENTREZID,
                     fun="enrichKEGG",
                     pvalueCutoff=0.2,
                     qvalueCutoff=0.5, 
                     minGSSize=0,
                     maxGSSize=1000)
ck <- ck%>%mutate(Cluster1=Cluster)
opfn2 <- paste(outdir, prefix, ".enrichKEGG.2.rds", sep="")
write_rds(ck, opfn2)
} ### 2, End                


#######################
### 3. Show figures ###
#######################

###3.1, MCl+treatment
if (FALSE){

cat("3.", "Show figures dotplots ordered by celltype+treatment", "\n") 

prefix <- "3_phiNew"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"
       
CL <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep="_")
lab2 <- setNames(gsub("-","+",CL), CL)
                      
###(1), showing GO
fn <- paste(outdir, prefix, ".enrichGO.rds", sep="")
cg <- read_rds(fn)  
cg2 <- cg%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))%>%filter(p.adjust<0.1)                 
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=15),
              axis.text.y=element_text(size=15))
        
figfn <- paste(outdir, prefix, ".Figure1.1_GO.png", sep="")
png(figfn, width=2500, height=2000, res=180)
print(fig1)
dev.off()
       
###(2), showing KEGG
fn <-  paste(outdir, prefix, ".enrichKEGG.rds", sep="")
ck <- read_rds(fn)
ck2 <- ck%>%mutate(ClusterNew=paste(contrast, MCls, sep="_"))%>%filter(p.adjust<0.1)     
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=12))
        
figfn <- paste(outdir, prefix, ".Figure1.2_KEGG.png", sep="")
png(figfn, width=1500, height=1500, res=150)
print(fig2)
dev.off()       
}

###3.2, cell type+treatment+direction
if (TRUE){

cat("3.2", "Show figures dotplots ordered by MCls+treat+direction", "\n")

prefix <- "3_phiNew"
outdir <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/"

cl <- paste( rep(c("LPS", "PHA", "LPS-DEX", "PHA-DEX"),each=4), 
            rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")
ii <- paste(rep(c("A","B","C","D"),each=4), rep(1:4,times=4), sep="")
cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".") 
cluster2 <- setNames(cl2, cl) 
lab2 <- setNames(gsub("-","+",cl),cl2)
                      
###(1)
fn <- paste(outdir, prefix, ".enrichGO.2.rds", sep="")
cg <- read_rds(fn) 
cg0 <- as.data.frame(cg)
x <- cluster2[as.character(cg0$Cluster)] 
cg2 <- cg%>%mutate(ClusterNew=x)%>%filter(p.adjust<0.1)                
fig1 <- enrichplot::dotplot(cg2, x=~ClusterNew, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=18),
              axis.text.y=element_text(size=12))
        
figfn <- paste(outdir, prefix, ".Figure2.1_GO.png", sep="")
png(figfn, width=4000, height=2600, res=180)
print(fig1)
dev.off()
       
###(2)
fn <-  paste(outdir, prefix, ".enrichKEGG.2.rds", sep="") 
ck <- read_rds(fn)
ck0 <- as.data.frame(ck)
x <- cluster2[as.character(ck0$Cluster)]  
ck2 <- ck%>%mutate(ClusterNew=x)%>%filter(p.adjust<0.1)     
fig2 <- dotplot(ck2, x=~ClusterNew, showCategory=10)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1,size=10),
              axis.text.y=element_text(size=12))
        
figfn <- paste(outdir, prefix, ".Figure2.2_KEGG.png", sep="")
png(figfn, width=1500, height=1500, res=150)
print(fig2)
dev.off()
       
}



