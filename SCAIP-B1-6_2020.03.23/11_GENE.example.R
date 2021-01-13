#
#rm(list=ls())
source("./Bin/LibraryPackage.R")

#outdir <- "./11_GENE.Example_output/tmp3_IFN/"
#if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)


########################
### defined function ###
########################

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
cg2 <- cg0%>%filter(grepl("type I interferon signaling pathway", Description))
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
df1 <- read_rds(fn)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))#%>%filter(abs(beta)>0.5, qval<0.1)
### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)%>%mutate(rn=paste(MCls, contrast, gene, sep="_"))

}###


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

#genels <- gene3[geneSel, "ENSEMBL"]
#symbolls <- gene3[geneSel, "SYMBOL"]
if(FALSE){
fn <- "./6_DEG.CelltypeNew_output/tmp1/2_meta.rds"
res <- read_rds(fn)%>%drop_na(beta,qval)#%>%filter(abs(beta)>0.5, qval<0.1)
gene2 <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

symbol_ls <- c("JAK1", "MX1", "MX2", "STAT1", "STAT2",
               "IRF1", "IRF7", "IRF8", "USP18", "OAS1", "OAS2", "OAS3", "OASL","IRF9", "TYK2", "IFNAR2", "IFNAR1")#, 
                #"FADD", "PTPN11", "METTL3")



lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
i <- 17          
symbol <- symbol_ls[i]
ens <- gene2%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
#ens <- ens[2]
oneMCl <- "Tcell"  

### expression
### fig 1
load("./6_DEG.CelltypeNew_output/tmp1/YtX.comb.RData")
load("./6_DEG.CelltypeNew_output/tmp1/0_ncell.RData")
#ddx <- dd%>%filter(ncell>20)
#YtX <- YtX[,ddx$bti]

count <- colSums(YtX)
count_after <- median(count)
count <- count/count_after
bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                  sampleID=cvt0[,3], Batch2=cvt0[,4])
                  
              
cvt2$y <- YtX[grepl(ens, rownames(YtX)),]/count 
d1 <- cvt2%>%drop_na(y)%>%filter(y>0)%>%filter(MCls==oneMCl)

fig1 <- ggplot(d1, aes(x=treats, y=log2(y), fill=treats))+
        geom_boxplot()+
        ylab(bquote(log[2]~"(Expression)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" expressed in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### dispersion
### fig2. box plot for dispersion
load("./10_RNA.Variance_output/tmp7/1_RNA.PhxNew.RData")
bti2 <- colnames(PhxNew2)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                   sampleID=cvt0[,3], Batch2=cvt0[,4])
                   
                   
cvt2$y <- PhxNew2[grepl(ens,rownames(PhxNew2)),]
d2 <- cvt2%>%drop_na(y)%>%filter(MCls==oneMCl)
fig2 <- ggplot(d2, aes(x=treats, y=log2(y), fill=treats))+
        geom_boxplot()+ylab(bquote(log[2]~"(Variability)"))+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" variability in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              

figfn <- paste("./11_GENE.Example_output/tmp3_IFN/Figure", i, "_", symbol, ".", oneMCl, ".png", sep="") 
#png(figfn, width=500, height=400, res=120)
png(figfn, width=700, height=600, res=150)
print(plot_grid(fig1, fig2, ncol=2)) 
dev.off()

} ###


###################################################
###          2, Example of gene                 ###
###  (1). DEG but not DVG; (2). DVG but not DEG ###
###################################################

if(FALSE){
load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
BG <- gsub("\\.[0-9].*", "", rownames(YtX_sel))
anno <- bitr(BG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)

### DEG
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
df1 <- read_rds(fn)
df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
x <- df1%>%group_by(gene)%>%summarise(n2=n(), .groups="drop")
gene2 <- x%>%filter(n2==4)%>%dplyr::pull(gene)
### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
x1 <- unique(df2$gene)
df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
x2 <- unique(df2$gene)

##
geneSel <- gene2[!gene2%in%x2]
gene
}###


###
### (1) plots, DEG but not DVG
if(FALSE){
outdir <- "./11_GENE.Example_output/4_DEGnoDVG/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
anno <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

symbol_ls <- c("CXCL8", "CXCL10", "IL6", "OAS2", "OAS3", "OASL")


lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
i <- 5          
symbol <- symbol_ls[i]
ens <- anno%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
#ens <- ens[2]
oneMCl <- "Tcell"  


### fig 1    
### expression
load("./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData")
load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)
YtX <- YtX[,ddx$bti]

count <- colSums(YtX)
count_after <- median(count)
count <- count/count_after
bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                  sampleID=cvt0[,3], Batch2=cvt0[,4])
                               
cvt2$y <- YtX[grepl(ens, rownames(YtX)),]/count 
d1 <- cvt2%>%drop_na(y)%>%filter(y>0)%>%filter(MCls==oneMCl)


fig1 <- ggplot(d1, aes(x=factor(treats), y=log2(y), fill=treats))+
        geom_boxplot(outlier.shape=1)+
        geom_signif(comparisons=list(c("CTRL", "LPS"), c("LPS", "LPS-DEX"),
                                     c("CTRL", "PHA"), c("PHA", "PHA-DEX")),
                    annotations=c("***", "***", "***", "***"),
                    y_position=c(9, 9.5, 10, 9.5), vjust=0.5)+    ##1,c(6.5,7,7.5,7)
        ylab(bquote(log[2]~"(Expression)"))+ylim(-1,10)+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" expressed in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### fig2. box plot for dispersion
load("./10_RNA.Variance_output/tmp9/1_RNA.PhxNew.RData")
bti2 <- colnames(PhxNew2)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                   sampleID=cvt0[,3], Batch2=cvt0[,4])
                                     
cvt2$y <- PhxNew2[grepl(ens,rownames(PhxNew2)),]
d2 <- cvt2%>%drop_na(y)%>%filter(MCls==oneMCl)
fig2 <- ggplot(d2, aes(x=treats, y=log2(y), fill=treats))+
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
} 

###
###

#load("./6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
#BG <- gsub("\\.[0-9].*", "", rownames(YtX_sel))
#anno <- bitr(BG, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb=org.Hs.eg.db)
#
#### DEG
#fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
#df1 <- read_rds(fn)
#df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
#x1 <- unique(df1$gene)
#### DVG
#fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
#df2 <- read.table(fn,header=T)
#df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
#x <- df2%>%group_by(gene)%>%summarise(n2=n(), .groups="drop")
#gene2 <- x%>%filter(n2==2)%>%dplyr::pull(gene)
#geneSel <- gene2[!gene2%in%x1]
#x2 <- anno%>%filter(ENSEMBL%in%geneSel)

### (2) plots, DVG but not DEG
if(FALSE){
outdir <- "./11_GENE.Example_output/5_DVGnoDEG/"
if (!file.exists(outdir)) dir.create(outdir, recursive=T, showWarnings=F)

fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
anno <- bitr(res$gene, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

symbol_ls <- c("RPL41", "RPL23A", "RPS27", "SUCO") 


lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
i <- 4         
symbol <- symbol_ls[i]
ens <- anno%>%filter(SYMBOL==symbol)%>%dplyr::pull(ENSEMBL)
#ens <- ens[2]
oneMCl <- "Tcell"  


### fig 1    
### expression
load("./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData")
load("./6_DEG.CelltypeNew_output/Filter2/0_ncell.RData")
ddx <- dd%>%filter(ncell>20)
YtX <- YtX[,ddx$bti]

count <- colSums(YtX)
count_after <- median(count)
count <- count/count_after
bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                  sampleID=cvt0[,3], Batch2=cvt0[,4])
                               
cvt2$y <- YtX[grepl(ens, rownames(YtX)),]/count 
d1 <- cvt2%>%drop_na(y)%>%filter(y>0)%>%filter(MCls==oneMCl)


fig1 <- ggplot(d1, aes(x=factor(treats), y=log2(y), fill=treats))+
        geom_boxplot(outlier.shape=1)+
        ylab(bquote(log[2]~"(Expression)"))+ylim(-1,5)+
        scale_fill_manual("", values=col1, labels=lab1)+
        scale_x_discrete("",labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" expressed in "~.(oneMCl)))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=-90, hjust=0, size=8),
              axis.title.y=element_text(size=10),
              plot.title=element_text(hjust=0.5, size=8),
              legend.position="none")
              
        
### fig2. box plot for dispersion
load("./10_RNA.Variance_output/tmp9/1_RNA.PhxNew.RData")
bti2 <- colnames(PhxNew2)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt2 <- data.frame(bti=bti2, MCls=cvt0[,1], treats=gsub("-EtOH", "", cvt0[,2]),
                   sampleID=cvt0[,3], Batch2=cvt0[,4])
                                     
cvt2$y <- PhxNew2[grepl(ens,rownames(PhxNew2)),]
d2 <- cvt2%>%drop_na(y)%>%filter(MCls==oneMCl)%>%mutate(rn2=paste(Batch2, treats, sep="_"))
tmp <- d2%>%group_by(rn2)%>%summarise(y=n(),.groups="drop")%>%filter(y>=3)
d2x <- d2%>%filter(rn2%in%tmp$rn2)


fig2 <- ggplot(d2x, aes(x=treats, y=log2(y), fill=treats))+
        geom_boxplot(outlier.shape=1)+ylab(bquote(log[2]~"(Variability)"))+
        geom_signif(comparisons=list(c("CTRL", "LPS"), #c("LPS", "LPS-DEX"),
                                     c("CTRL", "PHA")), #c("PHA", "PHA-DEX")),
                    annotations=c("***", "**"),
                    y_position=c(5, 6), tip_length = 0.01, vjust=0.5)+     ## 1,2 c(2,2.5,3,2.5)    
        scale_fill_manual("", values=col1, labels=lab1)+ylim(-6,6)+
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
} 




#################################################################
### (2). plots, spliced and unspliced reads along pseudotime  ###
#################################################################
           
if(FALSE){

lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

#fn <- "./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds"
#sc <- read_rds(fn)
    
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
for (oneMCl in MCls){

cat(oneMCl, "\n")


###
fn <- paste("./9_RNA.dynamic2_output/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="") 
meta <- read_rds(fn)
cellSel <- meta$NEW_BARCODE
sc2 <- subset(sc, cells=cellSel)

###
X <- sc2@assays$RNA@data
rnz <- rowSums(X)
X <- X[rnz>20,] 
rn <- rownames(X)

layers <- c("S", "U")
counts <- map(layers, function(layer){
   layer <- paste(layer, "-", sep="")
   rn1 <- rn[grepl(layer, rn)]
   x <- X[rn1,]
   counts <- colSums(x)
   counts_after <- median(counts)
   counts_after <- ifelse(counts_after==0, counts_after+1, counts_after)
   counts <- ifelse(counts!=0, counts/counts_after, counts+1)
})
names(counts) <- c("S", "U")

###loop by gene
for (i in 2:2){
gene0 <- genels[i]
symbol <- symbolls[i]
ii <- i 

dx <- try(map_dfr(layers, function(layer){
   gene <- paste(layer, gene0, sep="-")
   y0 <- X[grepl(gene,rn),]/counts[[layer]]   
   tmp <- data.frame(sampleID=meta$BEST.GUESS, treats=meta$treats, 
                     LDA_1=meta$LDA_1, 
                     MCls=meta$MCls, y=y0, index=layer)
}), silent=T)

if (class(dx)=="try-error") next

dx$treat2 <- gsub("-EtOH", "", dx$treats)

labUS <- c("S"="spliced", "U"="unspliced")

fig1 <- ggplot(dx%>%filter(y>0,index=="S"),aes(x=LDA_1, y=y, colour=treat2))+
        ylab("Normalized expression")+ 
        #ylab(bquote(log[2]~"(Expression)"))+
        geom_smooth(method="loess", se=F)+
        #scale_linetype_manual("", values=c("S"="solid", "U"="dotted"), labels=labUS)+
        scale_colour_manual("", values=col1, labels=lab1)+
        ggtitle(paste(symbol, " in ", oneMCl, sep=""))+
        theme_bw()+
        #facet_wrap(~MCls, scales="free_x")+
        theme(legend.background=element_rect(fill=NA),
              legend.text=element_text(size=8),
              plot.title=element_text(hjust=0.5, face="italic", size=8))
###              
figfn <- paste("./11_GENE.Example_output/tmp3_IFN/Figure", ii, ".3.1_", gene0, ".", oneMCl, ".fitting2.png", sep="")
png(figfn, width=550, height=400, res=120)
#png(figfn, width=500, height=400, res=130)
print(fig1) 
dev.off()
              
#fig2 <- ggplot(dx%>%filter(y>0,index=="S"), aes(x=LDA_1, y=log2(y)))+
#        #ylab("Normalized expression")+
#        ylab(bquote(log[2]~"(Expression)"))+
#        geom_point(aes(colour=treat2),size=0.2)+
#        scale_colour_manual("", values=col1, labels=lab1, guide=guide_legend(override.aes=list(size=1)))+
#        #facet_wrap(~index, nrow=2, labeller=labeller(index=labUS))+
#        theme_bw()+
#        theme(legend.background=element_blank(),
#              legend.text=element_text(size=8),
#              legend.key=element_rect(fill=NA),
#              legend.key.size=grid::unit(0.8,"lines"),
#              legend.position="none")
#fig2 <- ggMarginal(fig2, groupColour=T, groupFill=F, size=2)
              
#fig2 <- ggplot(dx%>%filter(y>0), aes(x=treat2, y=log2(y)))+
#        #ylab("Normalized expression")+
#        ylab(bquote(log[2]~"(Expression)"))+
#        geom_boxplot(aes(colour=treat2))+
#        scale_colour_manual("", values=col1, labels=lab1, guide=guide_legend(override.aes=list(size=1)))+
#        facet_wrap(~index, nrow=2, labeller=labeller(index=labUS))+
#        theme_bw()+
#        theme(legend.position="none")
                           
#figfn <- paste("./11_GENE.Example_output/tmp3_IFN/Figure", ii, ".3.2_", gene0, ".", oneMCl, ".S.png", sep="") 
#png(figfn, width=500, height=600, res=130)
#print(fig2) 
#dev.off()
   
} ###loop by gene

}###loop, MCls

} ###


#################
### (3),  Bin ###
#################
countFun <- function(X, layer){
   layer <- paste(layer, "-", sep="")
   rn <- rownames(X)
   rn1 <- rn[grepl(layer, rn)]
   x <- X[rn1,]
   counts <- colSums(x)
   counts_after <- median(counts)
   counts_after <- ifelse(counts_after==0, counts_after+1, counts_after)
   counts <- ifelse(counts!=0, counts/counts_after, counts+1)
   counts 
}

binFun <- function(x,breaks){
   L1<- cut(x, breaks,labels=1:breaks)
   d1 <- data.frame(x, L1)
}

plotData <- function(cvt, X, breaks=10, gene){
   cvt <- cvt%>%arrange(LDA)
   tmp <- binFun(cvt$LDA, breaks=breaks)       
   cvt <- cvt%>%mutate(Bin=tmp[,2])
   L1 <- tapply(cvt$LDA, cvt$Bin, mean)
   cvt$LDAave <- L1[cvt$Bin]
   
   Xi <- X[,cvt$NEW_BARCODE]
   rnz <- rowSums(Xi)
   #Xi <- Xi[rnz>20,]
   count_after <- countFun(Xi, "S")
   
   rn <- rownames(Xi)
   cvt$y <- Xi[grepl(paste("S-", gene, sep=""), rn),]     
   cvt
}


### Read data
if (FALSE){
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
fn <- "./5_IdenCelltype_output/4_SCAIP.MCls.Harmony.rds"
sc <- read_rds(fn)
X <- sc@assays$RNA@data


}

### plot along LDA_1
if(FALSE){
    
outdir <- "./11_GENE.Example_output/1_LDA1/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F)

fn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/7.1.LDA1_top100.rds"
topGene <- read_rds(fn)
oneMCl <- "Tcell"
top2 <- topGene%>%filter(MCls==oneMCl)
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="") 
meta <- read_rds(fn)

cvt <- meta%>%
       dplyr::select(NEW_BARCODE, BEST.GUESS, treats, LDA_1)%>%
       mutate(treat2=gsub("-EtOH", "", treats), MCls=oneMCl)%>%
       dplyr::rename(LDA=LDA_1)

i <- 3     
gene0 <- top2[i,"ENSEMBL"]
symbol <- top2[i,"SYMBOL"]

cvt2 <- plotData(cvt, X, breaks=10, gene=gene0) 

### (1). plot                
fig1 <- ggplot(cvt2%>%filter(y>0), aes(x=as.numeric(LDAave), y=y, group=Bin))+
        geom_boxplot()+
        #ylab("expression")+xlab("LDA_1")+
        #ylab(bquote(~log[2]~"(Normalized expression)"))+xlab("LDA_1")+ 
        ylab("Normalized data")+xlab("LDA_1")+
        ggtitle(bquote(~italic(.(symbol))~" in "~ .(oneMCl)))+
        theme_bw()+
        facet_wrap(~treat2, nrow=3, scales="free_x")+
        theme(plot.title=element_text(hjust=0.5, size=10),
              legend.position="none")
###              
figfn <- paste(outdir, "Figure4.", oneMCl,"_", i, ".", symbol, ".boxplot.png", sep="")
#png(figfn, width=550, height=400, res=120)
png(figfn, width=600, height=800, res=150)
print(fig1) 
dev.off()

###
### (2) plot curve
fig2 <- ggplot(cvt2%>%filter(y>0), aes(x=LDA, y=y, colour=treat2))+
        geom_smooth(method="loess", se=F)+
        #ylab("expression")+xlab("LDA_1")+
        #ylab(bquote(~log[2]~"(expression)"))+xlab("LDA_1")+
        ylab("Normalized data")+xlab("LDA_1")+
        scale_colour_manual("", values=col1, labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" in "~.(oneMCl)))+
        theme_bw()+
        theme(plot.title=element_text(hjust=0.5, size=10),
              legend.position="none")
              #legend.position=c(0.2,0.2))              
figfn <- paste(outdir, "Figure4.", oneMCl,"_", i, ".", symbol, ".fitting.png", sep="")
png(figfn, width=400, height=500, res=120)
#png(figfn, width=500, height=500, res=120)
print(fig2) 
dev.off()   
}


###LDA_2
if(TRUE){
    
outdir <- "./11_GENE.Example_output/2_LDA2/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
fn <- "./9_RNA.dynamic2_output/Filter2_DEG6571/Old/7.2.LDA2_top100.rds"
topGene <- read_rds(fn)
oneMCl <- "Tcell"
top2 <- topGene%>%filter(MCls==oneMCl)

###
fn <- paste("./9_RNA.dynamic2_output/Filter2_DEG6571/Old/4.1_LDA.", oneMCl, ".meta.rds", sep="") 
meta <- read_rds(fn)
                      
cvt <- meta%>%
       dplyr::select(NEW_BARCODE, BEST.GUESS, treats, LDA_2rev)%>%
       mutate(treat2=gsub("-EtOH", "", treats), MCls=oneMCl)%>%
       dplyr::rename(LDA=LDA_2rev)

i <- 3     
gene0 <- top2[i,"ENSEMBL"]
symbol <- top2[i,"SYMBOL"]

cvt2 <- plotData(cvt, X, breaks=10, gene=gene0) 

### (1). plot                
fig1 <- ggplot(cvt2%>%filter(y>0), aes(x=as.numeric(LDAave), y=y, group=Bin))+
        geom_boxplot()+
        #ylab("expression")+xlab("LDA_1")+
        #ylab(bquote(~log[2]~"(Normalized expression)"))+xlab("LDA_1")+ 
        ylab("Normalized data")+xlab("LDA_2")+
        ggtitle(bquote(~italic(.(symbol))~" in "~ .(oneMCl)))+
        theme_bw()+
        facet_wrap(~treat2, nrow=3, scales="free_x")+
        theme(plot.title=element_text(hjust=0.5, size=10),
              legend.position="none")
###              
figfn <- paste(outdir, "Figure4.", oneMCl,"_", i, ".", symbol, ".boxplot.png", sep="")
#png(figfn, width=550, height=400, res=120)
png(figfn, width=600, height=800, res=150)
print(fig1) 
dev.off()

###
### (2) plot curve
fig2 <- ggplot(cvt2%>%filter(y>0), aes(x=LDA, y=y, colour=treat2))+
        geom_smooth(method="loess", se=F)+
        #ylab("expression")+xlab("LDA_1")+
        #ylab(bquote(~log[2]~"(expression)"))+xlab("LDA_1")+
        ylab("Normalized data")+xlab("LDA_2")+
        scale_colour_manual("", values=col1, labels=lab1)+
        ggtitle(bquote(~italic(.(symbol))~" in "~.(oneMCl)))+
        theme_bw()+
        theme(plot.title=element_text(hjust=0.5, size=10),
              legend.position="none")
              #legend.position=c(0.2,0.2))              
figfn <- paste(outdir, "Figure4.", oneMCl,"_", i, ".", symbol, ".fitting.png", sep="")
png(figfn, width=400, height=500, res=120)
#png(figfn, width=500, height=500, res=120)
print(fig2) 
dev.off()   
}



