#
rm(list=ls())
source("./Bin/LibraryPackage.R")

#outdir <- "./11_GENE.Example_output/tmp3_IFN/"
#if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

###########################################
### Example genes show gene expression, ###
###  2020.11.4 last modified, by Julong ###
###########################################
if(FALSE){

#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#geneBG <- gsub("\\.[0-9].*", "", rownames(YtX))
#Bg <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
#cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
cg3 <- cg0%>%filter(grepl("type I interferon", Description))
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DVG
fn <- "./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew.enrichGO.rds"
cg <- read_rds(fn) 
cg0 <- as.data.frame(cg)
#cg3 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
#cg3 <- cg0%>%filter(grepl("regulation of cytokine production", Description))
#cg3 <- cg0%>%filter(Description=="cytokine receptor binding") 
cg3 <- cg0%>%filter(grepl("type I interferon", Description))
geneList <- unique(unlist(str_split(cg3$geneID,"/")))
gene3 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DEG
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
df1 <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))%>%filter(abs(beta)>0.5, qval<0.1)
### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))

}###

########################
### defined function ###
########################

###
adjGene <- function(cvt,center=T){
   cvt <- cvt%>%mutate(comb=paste(MCls, Batch, sep="_"))
   if(center){
      cvt <- cvt%>%group_by(comb)%>%mutate(yscale=y-mean(y,na.rm=T))%>%ungroup()
   }else{
      cvt <- cvt%>%mutate(yscale=y)
   }
}
###
adj2Gene <- function(cvt){
   cvt0 <- cvt%>%filter(treats=="CTRL")
   y0 <- cvt0$yscale
   names(y0) <- cvt0$sampleID
   y0[is.na(y0)] <- 0     
   cvt$y0 <- y0[cvt$sampleID]
   cvt <- cvt%>%mutate(yscale2=yscale-y0)
   cvt
}


###bulk, NB.mu, NB.phi
getData <- function(gene, datatype="bulk"){

### bulk data
if (datatype=="bulk"){
   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
   load(fn)
   rn <- gsub("\\..*", "", rownames(YtX_sel))
   rownames(YtX_sel) <- rn
   X <- YtX_sel[rn%in%gene,]+1

   bti <- colnames(YtX_sel)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])

   fn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
   load(fn)
   counts <- colSums(YtX)
   counts <- counts[colnames(YtX_sel)]

   X <- (X/counts)*1e+06
   cvt$y <- log2(X)
   cvt <- adjGene(cvt, center=T)
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
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
   cvt$y <- log2(X)
   cvt <- adjGene(cvt, center=T) 
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
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
   X <- log2(X)
   cvt$y <- X
   cvt <- adjGene(cvt, center=T) 
}
cvt
}
###
#######################################################################
### (1). boxplots show expression and dispersion for examples genes ###
#######################################################################

#genels1 <- c("ENSG00000172216", "ENSG00000152518", "ENSG00000135046", "ENSG00000120129",
#            "ENSG00000136689", "ENSG00000003402", "ENSG00000136244", "ENSG00000168209",
#            "ENSG00000128016", "ENSG00000170776", "ENSG00000073756", "ENSG00000157557",
#            "ENSG00000232810", "ENSG00000185650", "ENSG00000170345", "ENSG00000170266",
#            "ENSG00000138592", "ENSG00000119508", "ENSG00000153234", "ENSG00000187840")

#symbolls <- c("CEBPB", "ZFP36L2", "ANXA1", "DUSP1", "IL1RN", "CFLAR", "IL6",
#              "DDIT4", "ZFP36", "AKAP13", "PTGS2", "ETS2", "TNF", "ZFP36L1",
#               "FOS", "GLB1", "USP8", "NR4A3", "NR4A2", "EIF4EBP1") 

if(FALSE){

outdir <- "./11_GENE.Example_output/3_IFN/"
if ( !file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)
 
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%drop_na(beta,qval)#%>%filter(abs(beta)>0.5, qval<0.1)
gene2 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

#symbol_ls <- c("JAK1", "MX1", "MX2", "STAT1", "STAT2",
#               "IRF1", "IRF7", "IRF8", "USP18", "OAS1", "OAS2", "OAS3", "OASL","IRF9", "TYK2", "IFNAR2", "IFNAR1")#, 
#                #"FADD", "PTPN11", "METTL3")

lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
      
symbol_ls <- c("STAT1", "STAT2", "IRF9")           

i <- 3         
symbol <- symbol_ls[i]
ens <- gene2%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
oneMCl <- "Tcell"

cvt <- getData(gene=ens, datatype="NB.mu")%>%drop_na(yscale)%>%filter(MCls==oneMCl)
fig1 <- ggplot(cvt, aes(x=treats, y=yscale, fill=treats))+
        geom_boxplot()+
        ylab(bquote(log[2]~"(gene expression)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("", labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~"in T cell"))+
        theme_bw()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.x=element_text(angle=-90, hjust=0),
              legend.position="none")
figfn <- paste("./11_GENE.Example_output/3_IFN/Figure", i, "_", oneMCl, ".", symbol, ".png", sep="") 
png(figfn, width=300, height=400, res=120)
print(fig1) 
dev.off()             
        
### dispersion
### fig2. box plot for dispersion
#load("./10_RNA.Variance_output/tmp10/1_RNA.PhxNew.RData")
#bti2 <- colnames(PhxNew2)
#cvt0 <- str_split(bti2, "_", simplify=T)
#cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
#                   sampleID=cvt0[,3], Batch2=cvt0[,4])
                   
                   
#cvt2$y <- PhxNew2[grepl(ens,rownames(PhxNew2)),]
#d2 <- cvt2%>%drop_na(y)%>%filter(MCls==oneMCl)
#fig2 <- ggplot(d2, aes(x=treats, y=log2(y), fill=treats))+
#        geom_boxplot()+ylab(bquote(log[2]~"(Variability)"))+
#        scale_fill_manual("", values=col1, labels=lab1)+
#        scale_x_discrete("",labels=lab1)+
#        ggtitle(bquote(~italic(.(symbol))~" variability in "~.(oneMCl)))+
#        theme_bw()+
#        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
#              axis.title.y=element_text(size=10),
#              plot.title=element_text(hjust=0.5, size=8),
#              legend.position="none")
#              
#
#figfn <- paste("./11_GENE.Example_output/3_IFN/Figure", i, "_", symbol, ".", oneMCl, ".png", sep="") 
##png(figfn, width=500, height=400, res=120)
#png(figfn, width=700, height=600, res=150)
#print(plot_grid(fig1, fig2, ncol=2)) 
#dev.off()

} ###


###################################################
###          4, Example of gene                 ###
###           DMG but not DVG                   ###
###################################################

if(FALSE){     
###
outdir <- "./11_GENE.Example_output/4_DMGnoDVG/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
anno <- bitr(unique(res$gene), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df1 <- read.table(fn, header=T)
df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
x <- df1%>%group_by(gene)%>%summarise(n2=n(), .groups="drop")
gene1 <- x%>%filter(n2==4)%>%dplyr::pull(gene)
### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
gene2 <- unique(df2$gene)

##
geneSel <- gene1[!gene1%in%gene2]

anno2 <- anno%>%filter(ENSEMBL%in%geneSel)


lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
for (i in 1:nrow(anno2)){ 
       
symbol <- anno2$SYMBOL[i]
ens <- anno2$ENSEMBL[i]
oneMCl <- "Tcell"  
cat(i, symbol,"\n")
### fig 1, gene mean
cvt <- getData(gene=ens, datatype="NB.mu")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")    
fig1 <- ggplot(cvt2, aes(x=factor(treats), y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+
#        geom_signif(comparisons=list(c("CTRL", "LPS"), c("LPS", "LPS-DEX"),
#                                     c("CTRL", "PHA"), c("PHA", "PHA-DEX")),
#                    annotations=c("***", "***", "***", "***"),
#                    y_position=c(9, 9.5, 10, 9.5), vjust=0.5)+    ##1,c(6.5,7,7.5,7)
        ylab(bquote(log[2]~"(mean)"))+#ylim(-1,10)+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" mean in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### fig2. box plot for dispersion
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")
##
fig2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+ylab(bquote(log[2]~"(Variability)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" variability in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              

figfn <- paste(outdir, "Figure", i, "_", symbol, ".", oneMCl, ".png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=700, height=600, res=150)
print(plot_grid(fig1, fig2, ncol=2)) 
dev.off()
} ## loop by gene

} 


################################
###    5, Examples of genes  ###
###       DVG but not DMG    ###
################################

if(FALSE){     
###
outdir <- "./11_GENE.Example_output/5_DVGnoDMG/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
anno <- bitr(unique(res$gene), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df1a <- read.table(fn, header=T)
df1 <- df1a%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
gene1 <- unique(df1$gene)
### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
x <- df2%>%group_by(gene)%>%summarise(n2=n(),.groups="drop")
gene2 <- x%>%filter(n2==1)%>%dplyr::pull(gene)

##
geneSel <- gene2[!gene2%in%gene1]

anno2 <- anno%>%filter(ENSEMBL%in%geneSel)

lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###arrange effect size by decreasing order           
df3 <- res%>%
       filter(gene%in%geneSel,MCls=="Tcell")%>%
       arrange(desc(abs(beta)))%>%left_join(anno2, by=c("gene"="ENSEMBL"))  
symbol_ls <- unique(df3$SYMBOL)         
for (i in 1:20){ 
       
symbol <- symbol_ls[i]
ens <- anno2%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
oneMCl <- "Tcell"  
cat(i, symbol,"\n")
### fig 1, gene mean
cvt <- getData(gene=ens, datatype="NB.mu")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")    
fig1 <- ggplot(cvt2, aes(x=factor(treats), y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+
#        geom_signif(comparisons=list(c("CTRL", "LPS"), c("LPS", "LPS-DEX"),
#                                     c("CTRL", "PHA"), c("PHA", "PHA-DEX")),
#                    annotations=c("***", "***", "***", "***"),
#                    y_position=c(9, 9.5, 10, 9.5), vjust=0.5)+    ##1,c(6.5,7,7.5,7)
        ylab(bquote(log[2]~"(mean)"))+#ylim(-1,10)+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        scale_y_continuous(limits=c(-10,10),expand=expansion(mult=0.05))+
        ggtitle(bquote(~italic(.(symbol))~" mean in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### fig2. box plot for dispersion
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")
##
fig2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+ylab(bquote(log[2]~"(Variability)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" variability in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              

figfn <- paste(outdir, "Figure.gr1.", i, "_", symbol, ".", oneMCl, ".png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=700, height=600, res=150)
print(plot_grid(fig1, fig2, ncol=2)) 
dev.off()
} ## loop by gene

}


################################
### 5.1, shared DEG and DVG  ###
################################

if(TRUE){     
###
outdir <- "./11_GENE.Example_output/5.1_DVG.DMG/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
anno <- bitr(unique(res$gene), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df1a <- read.table(fn, header=T)
df1 <- df1a%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Monocyte")
x <- df1%>%group_by(gene)%>%summarise(n2=n(),.groups="drop")
gene1 <- x%>%filter(n2==4)%>%dplyr::pull(gene)

### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Monocyte")
x <- df2%>%group_by(gene)%>%summarise(n2=n(),.groups="drop")
gene2 <- x%>%filter(n2==4)%>%dplyr::pull(gene)

##
geneSel <- gene1[gene1%in%gene2]

anno2 <- anno%>%filter(ENSEMBL%in%geneSel)

lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###arrange effect size by decreasing order           
df3 <- df1a%>%
       filter(gene%in%geneSel, MCls=="Monocyte")%>%
       arrange(desc(abs(beta)))%>%left_join(anno2, by=c("gene"="ENSEMBL"))  
symbol_ls <- unique(df3$SYMBOL)         
for (i in 1:length(symbol_ls)){ 

i <- 6       
symbol <- symbol_ls[i]
ens <- anno2%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
oneMCl <- "Monocyte"  
cat(i, symbol,"\n")

### fig 1, gene mean                    NB.mu
cvt <- getData(gene=ens, datatype="NB.mu")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")    
fig1 <- ggplot(cvt2, aes(x=factor(treats), y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+
#        geom_signif(comparisons=list(c("CTRL", "LPS"), c("LPS", "LPS-DEX"),
#                                     c("CTRL", "PHA"), c("PHA", "PHA-DEX")),
#                    annotations=c("***", "***", "***", "***"),
#                    y_position=c(9, 9.5, 10, 9.5), vjust=0.5)+    ##1,c(6.5,7,7.5,7)
        ylab(bquote(log[2]~"(mean)"))+#ylim(-1,10)+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        #scale_y_continuous(limits=c(-10,10),expand=expansion(mult=0.05))+
        ggtitle(bquote(~italic(.(symbol))~" mean in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### fig2. box plot for dispersion
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")
##
fig2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
        geom_boxplot(outlier.shape=1)+ylab(bquote(log[2]~"(Variability)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" variability in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              

figfn <- paste(outdir, "Figure.gr4.", i, "_", symbol, ".", oneMCl, ".png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=700, height=600, res=150)
print(plot_grid(fig1, fig2, ncol=2)) 
dev.off()
} ## loop by gene

}


###############################################
#### 6, average expression in given pathway ###
###############################################
#load("./6_DEG.CelltypeNew_output/YtX.comb.RData")
#geneBG <- gsub("\\.[0-9].*", "", rownames(YtX))
#Bg <- bitr(geneBG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

##
avePathway <- function(X){
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


### 1
if (FALSE){
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
#cg2 <- cg0%>%filter(grepl("glucocorticoid", Description))
#cg2 <- cg0%>%filter(grepl("response to glucocorticoid", Description))
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
           
cvt <- getPathwayData(ens,datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
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
cvt <- getPathwayData(ens,datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
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
cvt <- getPathwayData(ens,datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
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
        
} ###



###(2)
if (FALSE){
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
#cg2 <- cg0%>%filter(grepl("myeloid leukocyte activation", Description))
cg2 <- cg0
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
           
cvt <- getPathwayData(ens,datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("myeloid leukocyte activation score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.1_myeloidleukocyteactivation.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens,datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("myeloid leukocyte activation score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.2_myeloidleukocyteactivation.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens,datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("myeloid leukocyte activation score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure2.3_myeloidleukocyteactivation.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        
} ###


###
### (3)
if (FALSE){
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("^response to glucocorticoid", Description))
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
           
cvt <- getPathwayData(ens,datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to glucocorticoid score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.1_glucocorticoid.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

###
cvt <- getPathwayData(ens,datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to glucocorticoid score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.2_glucocorticoid.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
cvt <- getPathwayData(ens,datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to glucocorticoid score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure3.3_glucocorticoid.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        
} ###


### 4
if (FALSE){
outdir <- "./11_GENE.Example_output/6_avePathway/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F) 

### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
#cg2 <- cg0%>%filter(grepl("^inflammatory response$", Description))
cg2 <- cg0%>%filter(grepl("^regulation of cytokine production$", Description))
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
cvt <- getPathwayData(ens,datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("regulation of cytokine production score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.1_cytokine.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()


### NB.mu
cvt <- getPathwayData(ens,datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("regulation of cytokine production score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.2_cytokine.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

### NB.phi
cvt <- getPathwayData(ens,datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("regulation of cytokine production score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure5.3_cytokine.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        
} ###


### 5
if (FALSE){
outdir <- "./11_GENE.Example_output/6_avePathway/"
### DEG
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("^response to bacterium$", Description))
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
cvt <- getPathwayData(ens,datatype="bulk")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig1 <- ggplot(cvt,aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.1_response.bacterium.bulk.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig1) 
dev.off()

### (2)
### NB.mu
cvt <- getPathwayData(ens,datatype="NB.mu")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig2 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.2_response.bacterium.NB.mu.png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=800, height=500, res=150)
print(fig2) 
dev.off()

###
### (3), NB.phi
cg <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew.enrichGO.2.rds") 
cg0 <- as.data.frame(cg)
cg2 <- cg0%>%filter(grepl("^response to bacterium$", Description))
geneList <- unique(unlist(str_split(cg2$geneID,"/")))
gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
ens <- gene2%>%dplyr::pull(ENSEMBL)
### NB.phi
cvt <- getPathwayData(ens,datatype="NB.phi")%>%mutate(treat2=gsub("-EtOH", "", treats))
fig3 <- ggplot(cvt%>%drop_na(y),aes(x=MCls, y=y, fill=treat2))+
        geom_boxplot()+
        ylab("response to bacterium score")+
        scale_fill_manual("", values=col1, labels=lab1)+
        theme_bw()+
        theme(axis.title.x=element_blank())
        
figfn <- paste(outdir, "Figure6.3_response.bacterium.NB.phi.png", sep="") 
png(figfn, width=800, height=500, res=150)
print(fig3) 
dev.off()
        
} ###







                                                                                   

 
