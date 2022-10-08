##
library(tidyverse)
library(parallel)
library(purrr)
library(reshape)
library(qqman)
library(qvalue)
##
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
library(ggsignif)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(viridis)
theme_set(theme_grey())

###
outdir <- "./10_RNA.Variance_output/tmp9_pub2/"
if ( !(file.exists(outdir))) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

 
######################
### used for paper ###
######################
rm(list=ls())
#### 1, barplots, DVG

Mybinom <- function(subdf){
   n1<- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
   n2 <- subdf%>%filter(direction==2)%>%dplyr::pull(ngene)
   ngene <- c(n1, n2)
   if(n1>n2){
      res <- binom.test(ngene, 0.5, alternative="greater")
   }else{
     res <- binom.test(ngene, 0.5, alternative="less") 
   }
   res$p.value
}
Mysymb <- function(pval){
  if(pval<0.001) symb <- "***"
  if(pval>=0.001 & pval<0.01) symb <- "**"
  if (pval>=0.01 & pval<0.05) symb <- "*"
  if(pval>0.05) symb <- ""
  symb
}
Mypos <- function(subdf){
 ny <- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
 ny
}

MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")

### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.5)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_") 
          
###read data
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)
x <- res%>%group_by(contrast)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))
x <- res%>%group_by(MCls)%>%nest()%>%mutate(ngene=map_dbl(data,~length(unique((.x)$gene))))

## up and down DGV
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2,-ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-400,300),5)
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
                          
## ###add star
## anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
##            mutate(pval=map_dbl(data, Mybinom), 
##                   symb=map_chr(pval, Mysymb),
##                   ypos=map_dbl(data, Mypos))%>%
##            unnest(cols=c(contrast,MCls))                          
                          
p0 <- ggplot(sig4, aes(x=MCls, y=ngene2))+
   geom_bar(aes(fill=comb),stat="identity")+
   scale_fill_manual(values=col2comb, labels="")+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
       vjust=ifelse(direction==2, 1.6, -0.5)), size=3.5)+ #
   scale_y_continuous(breaks=breaks_value, limits=c(-450,300),labels=abs(breaks_value))+
   ## ylab("Number of differentially variable genes")+    
   facet_grid(~contrast, labeller=facetlab)+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         ## axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12),
         strip.text.x=element_text(size=14),
         plot.margin=unit(c(5.5, 17, 5.5, 5.5), "points"))
              
## p1 <- p0+geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)        


figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.1_DVG.barplot_reviews.pdf"
## png(filename=figfn, width=850, height=400, res=120)
pdf(figfn, width=8.5, height=4)
print(p0)
grid.text("upregulated", x=unit(0.98,"npc"), y=unit(0.75,"npc"),
          rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
grid.text("downregulated", x=unit(0.98,"npc"), y=unit(0.4,"npc"),
          rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
dev.off()
## ggsave(figfn, width=8.5, height=4)


### binomial test between up and down regulated gene
## Mybinom <- function(subdf) {
##    n1<- subdf%>%filter(direction==1)%>%dplyr::pull(ngene)
##    n2 <- subdf%>%filter(direction==2)%>%dplyr::pull(ngene)
##    ngene <- c(n1, n2)
##    if(n1>n2){
##       res <- binom.test(ngene, 0.5, alternative="greater")
##    }else{
##      res <- binom.test(ngene, 0.5, alternative="less") 
##    }
##    res$p.value
## } 
## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T)%>%filter(qval<0.1, abs(beta)>0.5)

## ## up and down DGV
## sigs <- res%>%
##         mutate(direction=ifelse(beta>0, "1", "2"))%>%
##         group_by(contrast, MCls, direction)%>%
##         summarise(ngene=n(),.groups="drop")
## anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
##            mutate(pval=map_dbl(data,Mybinom))%>%
##            unnest(cols=c(contrast,MCls))
           



#####################
### heatmap plots ###
#####################
rm(list=ls())
load("./10_RNA.Variance_output/tmp9/Sig3x.DGP.RData")
Geneunq <- sigs

col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

######################
### (1) Using beta ###
######################  
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta2"
res <- read.table(fn,header=T)
TMP <- map_dfc(MCls, function(oneMCl){
   tmp <- map_dfc(Contrast, function(oneC){
      d0 <- res %>% filter(MCls==oneMCl, contrast==oneC)
      rownames(d0) <- d0$gene
      beta0 <- d0[Geneunq,"beta"]
      beta0
   })
   tmp 
})
rownames(TMP) <- Geneunq
TMP <- as.data.frame(TMP)
contrast2 <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
conditions <- paste(rep(MCls,each=4), rep(contrast2, times=4), sep="_")
names(TMP) <- conditions
ngene <- nrow(TMP)

###(1) heatmap
###mybreaks
ii <- rowSums(is.na(TMP))
TMP0 <- TMP[ii==0,]

### (2) correlation heatmap ###
#Neworder <- c("Bcell_LPS+DEX", "Bcell_LPS", "Bcell_PHA+DEX",  "Bcell_PHA",
#              "Monocyte_LPS+DEX", "Monocyte_LPS", "Monocyte_PHA+DEX", "Monocyte_PHA", 
#              "Tcell_LPS+DEX","Tcell_LPS", "Tcell_PHA+DEX", "Tcell_PHA", 
#               "NKcell_LPS+DEX", "NKcell_LPS", "NKcell_PHA+DEX", "NKcell_PHA") 
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "NKcell_LPS+DEX", "Tcell_LPS+DEX",
#              "NKcell_PHA+DEX", "Tcell_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",               
#              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#              "NKcell_PHA", "Tcell_PHA","NKcell_LPS", "Tcell_LPS") 
Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
              "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
              "Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "NKcell_LPS+DEX", "NKcell_PHA+DEX", "Tcell_LPS+DEX", "Tcell_PHA+DEX")

corr <- cor(TMP0, method="spearman")[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)

x <- str_split(colnames(corr), "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], contrast=x[,2])
rownames(tmp_column) <- colnames(corr)
tmp_colors <- list(celltype=col2, contrast=col1)

p2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
   ## cellwidth=18, cellheight=18,            
   cluster_rows=F, cluster_cols=F, annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T, annotation=, fontsize_row=12,
   show_colnames=T, show_rownames=F, na_col="white")

p2 <- ggplotify::as.ggplot(p2)

## write_rds(fig2, file="./10_RNA.Variance_output/tmp9_pub/3.3_corr.beta.rds")
figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.2_corr.beta.png"
png(figfn, width=600, height=600, res=120)
print(p2)
dev.off()


#######################################################
#### 3, scatter plot of beta between va, phi and mu ###
#######################################################
rm(list=ls())       
### my fun 1, generate data frame used for plots 
myDFxy <- function(dfx, dfy){
###
   dfx <- dfx%>%dplyr::select(zscore, beta, qval, rn, contrast, MCls, gene)
   dfy <- dfy%>%dplyr::select(zscore, beta, qval, rn)       
   dfxy <- dfx%>%inner_join(dfy, by="rn")
###        
   x <- dfxy$qval.x
   y <- dfxy$qval.y
   Bx <- abs(dfxy$beta.x)
   By <- abs(dfxy$beta.y)
   gr <- rep(1, nrow(dfxy))
   gr[(x<0.1&Bx>0.5)&((y>=0.1)|(y<0.1&By<=0.5))] <- 2
   gr[((x>=0.1)|(x<0.1&Bx<=0.5))&(y<0.1&By>0.5)] <- 3
   gr[(x<0.1&Bx>0.5)&(y<0.1&By>0.5)] <- 4
   dfxy$gr <- gr
###
   dfxy
}

### label expression function       
feq <- function(x){
  #r <- format(as.numeric(x$estimate),digits=1)
  r <- round(as.numeric(x$estimate), digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  eq
  #r 
}

#################
### Read data ###
#################

### df2, residual dispersion
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(file=fn, header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 
       
### df3, mean expression
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df3 <- read.table(file=fn, header=T)%>%drop_na(beta,qval)%>%mutate(zscore=beta/stderr) 

###
mycol <- c("1"="grey", "2"="#118C4F", "3"="blue", "4"="#9400D3")
#mycol <- c("1"="#bababa", "2"="#de2d26", "3"="#6baed6", "4"="#756bb1")
Newcon2 <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", "PHA"="PHA", "PHA-DEX"="PHA+DEX")


#### (1). mean vs residual dispersion
dfxy1 <- myDFxy(df3, df2)
mylabel <- c("1"="NS", "2"="DEG(only)", "3"="DVG(only)", "4"="Both")

anno_df1 <- dfxy1%>%
   group_by(contrast, MCls)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")), r=map_dbl(data,~cor((.x)$zscore.x, (.x)$zscore.y)),
   eq=map(corr,feq))%>%
   dplyr::select(-data,-corr)

###scatterplots
dfxy2 <- dfxy1%>%filter(MCls=="Monocyte",grepl("PHA-DEX", contrast))
anno2 <- anno_df1%>%filter(MCls=="Monocyte", contrast=="PHA-DEX")
p1 <- ggplot(dfxy2, aes(x=zscore.x, y=zscore.y))+
   geom_point(aes(colour=factor(gr)), size=0.3)+
   scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
   geom_hline(yintercept=0, color="grey60")+
   geom_vline(xintercept=0, color="grey60")+
   geom_text(data=anno2, x=32, y=-3, aes(label=eq), size=3.5, parse=T)+
   ggtitle("Monocytes PHA+DEX")+ 
   xlab("zscore of differential gene expression")+ylab("zscore of differential gene variability")+
      ## facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.4,units="cm"),
         legend.position=c(0.8, 0.8),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         axis.text=element_text(size=10),
         axis.title=element_text(size=14),
         plot.title=element_text(hjust=0.5,size=14))

###
###
## anno3 <- anno_df1%>%dplyr::select(-eq)%>%as.data.frame()
## anno3$x <- 1
## p2 <- ggplot(anno3,aes(x=x, y=r))+
##    geom_violin(width=0.5)+
##    geom_jitter(aes(color=factor(MCls)), width=0.08, size=0.25)+
##    scale_color_manual(values=c(c("Bcell"="#4daf4a",
##        "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")))+ 
##    ylim(-0.5,0.2)+
##    geom_hline(yintercept=0, color="grey60")+
##    ggtitle("rho")+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.text.x=element_blank(),
##          axis.text.y=element_text(size=8),
##          axis.title=element_blank(),
##          axis.ticks.x=element_blank(),
##          plot.title=element_text(size=10, hjust=0.5, vjust=-1))

## fig0 <- p1+
##    inset_element(p2, left=0.6, bottom=0.6,
##        right=unit(1,"npc")-unit(0.2,"cm"), top=unit(1,"npc")-unit(0.2,"cm"))


###              
figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.3_DMGvsDVG.png"
png(filename=figfn, width=500, height=420, res=120)  
print(p1)
dev.off()


###
###
anno_df1 <- dfxy1%>%
   group_by(MCls)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore.x, (.x)$zscore.y, method="pearson")), r=map_dbl(data,~cor((.x)$zscore.x, (.x)$zscore.y)),
   eq=map(corr,feq))%>%
   dplyr::select(-data,-corr)
###
dfxy2 <- dfxy1%>%filter(MCls=="Monocyte")
anno2 <- anno_df1%>%filter(MCls=="Monocyte")

p0 <- ggplot(dfxy2, aes(x=zscore.x, y=zscore.y))+
   geom_point(aes(colour=factor(gr)), size=0.3)+
   scale_color_manual(values=mycol, labels=mylabel, guide=guide_legend(override.aes=list(size=1.5)))+
   geom_hline(yintercept=0, color="grey60")+
   geom_vline(xintercept=0, color="grey60")+
   geom_text(data=anno2, x=40, y=-10, aes(label=eq), size=3.5, parse=T)+
   ggtitle("Monocyte")+ 
   xlab("zscore of differential gene expression")+ylab("zscore of differential gene variability")+
      ## facet_grid(contrast~MCls, labeller=labeller(contrast=Newcon2))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.text=element_text(size=9),
         legend.key.size=unit(0.4,units="cm"),
         legend.position=c(0.8, 0.8),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         axis.text=element_text(size=10),
         axis.title=element_text(size=12),
         plot.title=element_text(hjust=0.5,size=14))

figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.3.2_DMGvsDVG.png"
png(filename=figfn, width=500, height=420, res=120)  
print(p0)
dev.off()

###sumamary number
x <- dfxy1%>%
     filter(gr==4)%>%
     group_by(MCls, contrast)%>%nest()%>%
     mutate(ngene=map_dbl(data,nrow))



#####################
### example genes ###
#####################

rm(list=ls())

########################
### defined function ###
########################

adjGene <- function(cvt, center=T){
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


### 3.4, Example of gene, DMG but not DVG
####

## fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
## res <- read.table(fn, header=T)
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
geneAll <- unique(read_rds(fn)$gene)
anno <- bitr(geneAll, fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
## fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
## df1 <- read.table(fn, header=T)
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
df1 <- read_rds(fn)
df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Tcell")
pval <- df1$pval
symbol <- rep("ns", nrow(df1))
symbol[pval<=0.05] <- "*"
symbol[pval<=0.01] <- "**"
symbol[pval<=0.001] <- "***"
df1$symbol <- symbol
                  
###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
### 
anno2 <- anno%>%filter(SYMBOL=="TAP2")
symbol <- anno2[1,"SYMBOL"]
ens <- anno2[1, "ENSEMBL"]    
oneMCl <- "Tcell"  

### fig 1, gene mean
cvt <- getData(gene=ens, datatype="bulk")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")

### annotation
sig_df <- df1%>%filter(gene==ens)%>%dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2), .groups="drop")
sig_df <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast")  

p1 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+    
   geom_text(data=sig_df,
      aes(x=treats, y=ymax+0.2, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~.(oneMCl)~"("~italic(.(symbol))~")"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_blank(),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=14),
         legend.position="none")


### fig2. box plot for dispersion
## plot data
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")
###
### annotation data
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
df2 <- df2%>%filter(MCls==oneMCl, gene==ens)

##
p2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+    
   ylab(bquote(~log[2]~"(Variability)"))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete(labels=lab1)+
   ## ggtitle(bquote(~.(oneMCl)))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         ## plot.title=element_text(hjust=0.5, size=12),
         plot.title=element_blank(),
         legend.position="none")

## fig_ls <- list(fig1=fig1, fig2=fig2)
## write_rds(fig_ls, "./10_RNA.Variance_output/tmp9_pub/3.4_onlyDMG.rds")


figfn <- paste("./10_RNA.Variance_output/tmp9_pub/Figure3.4.2_",
    oneMCl,".",symbol, ".png", sep="")
png(figfn, width=480, height=620, res=120)
print(plot_grid(p1, p2, nrow=2,
   labels=NULL, label_size=16, label_x=0,             
   label_fontface="plain", align="v", axis="lr",             
   rel_heights=c(1,1.2))) 
dev.off()



####
###

## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## for (i in 1:4){
## ##    
## oneMCl <- MCls[i]
## geneName <- "CD52"
## cat(oneMCl, geneName, "\n")    

## fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
## df1 <- read_rds(fn)
## df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls==oneMCl)
## pval <- df1$pval
## symbol <- rep("ns", nrow(df1))
## symbol[pval<=0.05] <- "*"
## symbol[pval<=0.01] <- "**"
## symbol[pval<=0.001] <- "***"
## df1$symbol <- symbol
                  
## ###
## lab1 <- c("CTRL"="CTRL", 
##           "LPS"="LPS", "LPS-DEX"="LPS+DEX",
##           "PHA"="PHA", "PHA-DEX"="PHA+DEX")
## col1 <- c("CTRL"="#828282", 
##            "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##            "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
           
## ### 
## anno2 <- anno%>%filter(SYMBOL==geneName)
## symbol <- anno2[1,"SYMBOL"]
## ens <- anno2[1, "ENSEMBL"]    
## ## oneMCl <- "Tcell"  

## ### fig 1, gene mean
## cvt <- getData(gene=ens, datatype="bulk")%>%filter(MCls==oneMCl)
## cvt <- adj2Gene(cvt)
## cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")

## ### annotation
## sig_df <- df1%>%filter(gene==ens)%>%dplyr::select(MCls, contrast, symbol)
## pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2, na.rm=T), .groups="drop")
## sig_df <- sig_df%>%
##    left_join(pos_df, by=c("contrast"="treats"))%>%
##    dplyr::rename("treats"="contrast") 

## p1 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
##    geom_violin(width=0.8)+
##    geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+    
##    geom_text(data=sig_df,
##       aes(x=treats, y=ymax+0.2, label=symbol))+
##    ylab(bquote(~log[2]~"(Expression)"))+
##    scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
##    scale_fill_manual("", values=col1, labels=lab1)+
##    scale_x_discrete("", labels=lab1)+
##    ggtitle(bquote(~.(oneMCl)~"("~italic(.(symbol))~")"))+
##    theme_bw()+
##    theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
##          ## axis.text.x=element_blank(),
##          axis.title.x=element_blank(),
##          axis.ticks.x=element_blank(),
##          axis.title.y=element_text(size=12),
##          plot.title=element_text(hjust=0.5, size=14),
##          legend.position="none")
## figfn <- paste("./10_RNA.Variance_output/tmp9_pub2/Figure3.", i, "_", oneMCl, "_", symbol, ".png", sep="")
## png(figfn, width=420, height=420, res=120)
## print(p1) 
## dev.off()
## }###
## fig1 <- plot_grid(p1, p2, nrow=2,
##    align="v", axis="lr", rel_heights=c(1,1.2))


###
###  3.5, Examples of genes, DVG but not DMG

fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
anno <- bitr(unique(res$gene), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df1 <- read.table(fn, header=T)%>%filter(MCls=="Tcell")

### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn, header=T)%>%filter(MCls=="Tcell") 
pval <- df2$pval
symbol <- rep("ns", nrow(df2))
symbol[pval<=0.05] <- "*"
symbol[pval<=0.01] <- "**"
symbol[pval<=0.001] <- "***"
df2$symbol <- symbol


lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

oneMCl <- "Tcell"
anno2 <- anno%>%filter(SYMBOL=="RPL23A")       
ens <- anno2[1,"ENSEMBL"]
symbol <- anno2[1,"SYMBOL"]

### fig 1, gene mean
cvt <- getData(gene=ens, datatype="NB.mu")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")
sig_df1 <- df1%>%filter(gene==ens)

p1 <- ggplot(cvt2, aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+ 
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   ylab(bquote(~log[2]~"(expression)"))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete(labels=lab1)+ 
   scale_y_continuous(limits=c(-1,1))+
   ggtitle(bquote(~.(oneMCl)~"("~italic(.(symbol))~")"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=14),
         legend.position="none")
              
        
### fig2. box plot for dispersion
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")

## annotation
sig_df <- df2%>%filter(gene==ens)%>%dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2), .groups="drop")
sig_df <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast") 

p2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
   geom_violin(width=0.8)+ 
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_text(data=sig_df, aes(x=treats, y=ymax+0.2, label=symbol), size=3)+    
   ylab(bquote(~log[2]~"(Variability)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2,0.2)))+    
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("",labels=lab1)+
   ## ggtitle(bquote(~italic(.(symbol))~" in "~.(oneMCl)))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         ## plot.title=element_text(hjust=0.5, size=12),
         plot.title=element_blank(),
         legend.position="none")

## fig_ls <- list(fig1=fig1, fig2=fig2)
## write_rds(fig_ls, file="./10_RNA.Variance_output/tmp9_pub/3.5_onlyDVG.rds")
###
figfn <- paste("./10_RNA.Variance_output/tmp9_pub/Figure3.5.2_",
    oneMCl,".",symbol, ".png", sep="")
## png(figfn, width=620, height=450, res=120)
## print(plot_grid(p1, p2, ncol=2)) 
## dev.off()
png(figfn, width=480, height=620, res=120)
print(plot_grid(p1, p2, nrow=2,
   labels=NULL, label_size=16, label_x=0,             
   label_fontface="plain", align="v", axis="lr",             
   rel_heights=c(1,1.2))) 
dev.off()

## fig2 <- plot_grid(p1, p2, nrow=2,
##    align="v", axis="lr", rel_heights=c(1,1.2))


###
### 3.6, shared DEG and DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
res <- read.table(fn, header=T)
anno <- bitr(unique(res$gene), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL"), OrgDb=org.Hs.eg.db)

### DMG
fn <- "./10_RNA.Variance_output/tmp9/4_mu.meta"
df1 <- read.table(fn, header=T)%>%filter(MCls=="Monocyte")
df1 <- df1%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Monocyte")
###
pval <- df1$pval
symbol <- rep("ns", nrow(df1))
symbol[pval<=0.05] <- "*"
symbol[pval<=0.01] <- "**"
symbol[pval<=0.001] <- "***"
df1$symbol <- symbol

### DVG
fn <- "./10_RNA.Variance_output/tmp9/3_phiNew.meta"
df2 <- read.table(fn,header=T)
df2 <- df2%>%filter(abs(beta)>0.5,qval<0.1)%>%filter(MCls=="Monocyte")
##
pval <- df2$pval
symbol <- rep("ns", nrow(df2))
symbol[pval<=0.05] <- "*"
symbol[pval<=0.01] <- "**"
symbol[pval<=0.001] <- "***"
df2$symbol <- symbol

##
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###
oneMCl <- "Monocyte"
anno2 <- anno%>%filter(SYMBOL=="CCL3L1")
ens <- anno2[1, "ENSEMBL"]
symbol <- anno2[1,"SYMBOL"]  

### fig 1, gene mean, NB.mu
cvt <- getData(gene=ens, datatype="NB.mu")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y, yscale2)%>%filter(treats!="CTRL")

sig_df <- df1%>%filter(gene==ens)%>%dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2), .groups="drop")
sig_df1 <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast")

###
p1 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+ 
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_text(data=sig_df1,
             aes(x=treats, y=ymax+0.5, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete(labels=lab1)+
   ggtitle(bquote(~.(oneMCl)~"("~italic(.(symbol))~")"))+    
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         plot.title=element_text(hjust=0.5, size=14),
         legend.position="none")
    
                      
### fig2. box plot for dispersion
cvt <- getData(gene=ens, datatype="NB.phi")%>%filter(MCls==oneMCl)
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y,yscale2)%>%filter(treats!="CTRL")
    
###
sig_df <- df2%>%filter(gene==ens)%>%dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2), .groups="drop")
sig_df2 <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast")

p2 <- ggplot(cvt2, aes(x=treats, y=yscale2, fill=treats))+
   geom_violin(width=0.8)+ 
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_text(data=sig_df2,
      aes(x=treats, y=ymax+0.2, label=symbol), size=3)+    
   ylab(bquote(~log[2]~"(Variability)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2,0.2)))+    
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete(labels=lab1)+
   ## ggtitle(bquote(~.(oneMCl)))+
   theme_bw()+
   theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         ## plot.title=element_text(hjust=0.5, size=12),
         plot.title=element_blank(),
         legend.position="none")

## fig_ls <- list(fig1=fig1, fig2=fig2)
## write_rds(fig_ls, "./10_RNA.Variance_output/tmp9_pub/3.6_Both.rds")

figfn <- paste("./10_RNA.Variance_output/tmp9_pub/Figure3.6.2_",
    oneMCl,".",symbol, ".png", sep="")
## png(figfn, width=620, height=450, res=120)
png(figfn, width=480, height=620, res=120)
print(plot_grid(p1, p2, nrow=2,
   ## labels=c("F", ""), label_size=16, label_fontface="plain",
   align="v", axis="lr",             
   rel_heights=c(1,1.2))) 
dev.off()

## fig3 <- plot_grid(p1, p2, nrow=2,
##     align="v", axis="lr", rel_heights=c(1,1.2))


###
## figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.4_comb.png"
## png(figfn, width=1200, height=500, res=120)
## print(plot_grid(fig1, fig2, fig3, ncol=3,
##    labels=c("D", "E", "F"), label_x=c(0.05, 0.05, 0.1),
##    label_size=16, label_fontface="plain"))
## dev.off()



