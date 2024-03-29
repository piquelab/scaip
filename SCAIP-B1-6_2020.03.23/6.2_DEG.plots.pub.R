##
###
rm(list=ls())
###
library(tidyverse)
library(parallel)
##
library(DESeq2)
library(annotables)
library(biobroom)
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggplotify)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(ComplexHeatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(circlize)

## theme_set(theme_grey())


outdir <- "./6_DEG.CelltypeNew_output/Filter2_pub2/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


###
### barplots showing DEG across cell types and contrast 
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
 
### facet by MCls
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
### colors
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
col2w <- colorspace::lighten(col2,0.5)
col2comb <- c(col2,col2w)
names(col2comb) <- paste(MCls, rep(c(1,2),each=4), sep="_")

### read data
fn <- "./6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)%>%filter(qval<0.1,abs(beta)>0.5)%>%drop_na(beta)
sigs <- res%>%
        mutate(direction=ifelse(beta>0, "1", "2"))%>%
        group_by(contrast, MCls, direction)%>%
        summarise(ngene=n(),.groups="drop")
        
sig4 <- sigs%>%mutate(ngene2=ifelse(direction==2, -ngene, ngene),
                      comb=paste(MCls, direction, sep="_"))
breaks_value <- pretty(c(-1800,1800),5)
facetlab <- as_labeller(c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
                          "PHA"="PHA", "PHA-DEX"="PHA+DEX"))
                          
###add star
## anno_df <- sigs%>%group_by(contrast, MCls)%>%nest()%>%
##            mutate(pval=map_dbl(data, Mybinom), 
##                   symb=map_chr(pval, Mysymb),
##                   ypos=map_dbl(data, Mypos))%>%
##            unnest(cols=c(contrast,MCls))                          

fig0 <- ggplot(sig4, aes(x=MCls, y=ngene2))+
   geom_bar(aes(fill=comb),stat="identity")+
   scale_fill_manual(values=col2comb, labels="")+
   geom_hline(yintercept=0, color="grey60")+
   geom_text(aes(x=MCls, y=ngene2, label=abs(ngene2), 
      vjust=ifelse(direction==2, 1.5, -0.5)), size=3.5)+        
   scale_y_continuous(breaks=breaks_value,
      limits=c(-2000,2000),labels=abs(breaks_value))+
   ## ylab("Number of differentially expressed genes")+ 
   facet_grid(~contrast, labeller=facetlab)+   
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         ## axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=12),
         axis.text.y=element_text(size=12),
         strip.text.x=element_text(size=14),
         plot.margin=unit(c(5.5, 17, 5.5, 5.5), "points"))
## p1 <- fig0+
##    geom_text(data=anno_df, aes(x=MCls, y=ypos, label=symb), colour="black", vjust=-1, size=3)


## png("./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.1_barplot_reviews.png", width=850, height=400, res=120)
pdf("./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.1_barplot_reviews.pdf", width=8.5, height=4)
print(fig0)
grid.text("upregulated", x=unit(0.98,"npc"), y=unit(0.7,"npc"),
          rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
grid.text("downregulated", x=unit(0.98,"npc"), y=unit(0.38,"npc"),
          rot=90, hjust=0.5, vjust=0.5, gp=gpar(cex=0.9))
dev.off()

###############
### heatmap ###
###############

load("./6_DEG.CelltypeNew_output/Filter2/Sigs.gene.DEG.RData")
DEGunq <- unique(sigs)

col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
          "NKcell"="#aa4b56", "Tcell"="#ffaa00")

###
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
fn <- "./6_DEG.CelltypeNew_output/Filter2/3_Batch1456.meta.rds"
res <- read_rds(fn)
TMP <- map_dfc(MCls, function(oneMCl){
   tmp <- map_dfc(Contrast, function(oneC){ 
      d0 <- res %>% filter(MCls==oneMCl, contrast==oneC)%>%mutate(z=beta/stderr)
      rownames(d0) <- d0$gene
      beta0 <- d0[DEGunq,"beta"]
      beta0
   })
   tmp 
})
rownames(TMP) <- DEGunq
TMP <- as.data.frame(TMP)
contrast2 <- c("LPS", "LPS+DEX", "PHA", "PHA+DEX")
conditions <- paste(rep(MCls,each=4), rep(contrast2, times=4), sep="_")
names(TMP) <- conditions
ngene <- nrow(TMP)

ii <- rowSums(is.na(TMP))
TMP0 <- TMP[ii==0,]



###(1) heatmap
###mybreaks
y <- do.call(c, TMP0)
y0 <- y[abs(y)<2] #99% percent quantile(abs(y),probs=0.99)
mybreaks <- c(min(y),quantile(y0,probs=seq(0,1,length.out=98)),max(y))
names(mybreaks) <- NULL

###
### new order
Neworder <- c("Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA",
              "NKcell_LPS", "NKcell_PHA", "Tcell_LPS", "Tcell_PHA",
              "Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
              "NKcell_LPS+DEX", "NKcell_PHA+DEX", "Tcell_LPS+DEX", "Tcell_PHA+DEX")
TMP0 <- TMP0[,Neworder]


###colors pheatmap parameter
## x <- str_split(Neworder, "_", simplify=T)
## tmp_column <- data.frame("celltype"=x[,1], "contrast"=x[,2])
## rownames(tmp_column) <- Neworder
## tmp_colors <- list("celltype"=col2, contrast=col1) #brewer.pal(4,"Set1")

## mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)

TMP2 <- t(TMP0)
#mycol <- viridisLite::viridis(100)
#mycol <- viridisLite::cividis(100, direction=1)
## fig1 <- pheatmap(TMP2, col=mycol, breaks=mybreaks, 
##          scale="none",
##          border_color="NA",
##          cluster_rows=T, cluster_cols=T, 
##          annotation_row=tmp_column,
##          annotation_colors=tmp_colors, annotation_legend=F,
##          show_colnames=F, show_rownames=T, fontsize=12,
##          na_col="white")
 
## figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.2_1_heatmap.beta.png"
## png(figfn, width=1100, height=550,res=120)
## print(fig1)
## dev.off()
x <- str_split(Neworder, "_", simplify=T)
anno_df <- data.frame("contrast"=x[,2],"celltype"=x[,1])
 
row_ha <- rowAnnotation(df=anno_df,
    col=list(contrast=col1, celltype=col2),
    show_legend=c(celltype=F, contrast=F),
    annotation_name_side="top",
    show_annotation_name=T)



mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)

mycol2 <- colorRamp2(mybreaks, mycol)




fig1 <- Heatmap(TMP2, col=mycol2,
   cluster_rows=T, cluster_columns=T,
   show_row_names=T, row_names_side="right",
   row_names_gp=gpar(fontsize=10),
   show_row_dend=T, show_column_dend=T,
   show_column_names=F,
   left_annotation=row_ha,
   ##
   heatmap_legend_param=list(title="LFC",
      title_gp=gpar(fontsize=12),
      labels_gp=gpar(fontsize=10),
      grid_width=grid::unit(0.6, "cm"),
      legend_height=grid::unit(5, "cm")),
  ## 
   use_raster=T, raster_device="png")


figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.2_1_heatmap.beta.png"
png(figfn, width=1100, height=580,res=120)
set.seed(0)
fig1 <- draw(fig1)
r.list <- row_order(fig1)
r.dend <- row_dend(fig1)
dev.off()





### (2) correlation heatmap ###
#Neworder <- c("Monocyte_LPS+DEX", "Monocyte_PHA+DEX", "Bcell_LPS+DEX", "Bcell_PHA+DEX",
#              "Tcell_LPS+DEX", "Tcell_PHA+DEX", "NKcell_LPS+DEX", "NKcell_PHA+DEX",
#              "Monocyte_LPS", "Monocyte_PHA", "Bcell_LPS", "Bcell_PHA", 
#               "Tcell_LPS", "Tcell_PHA", "NKcell_LPS", "NKcell_PHA") 
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
   cellwidth=18, cellheight=18,
   cluster_rows=F, cluster_cols=F,
   annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
   show_colnames=T, show_rownames=F, na_col="white", fontsize=10, fontsize_row=12)
## p2 <- as.ggplot(p2)

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.2_corr.beta.png"
png(figfn, width=720, height=720,res=120)
print(p2)
dev.off()




########################
### enriched pathway ###
########################

odds.fun <-  function(df){
###    
   res <- map_dfr(1:nrow(df), function(i){
      Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
      Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
      dat <- data.frame(Diff=Diff, Bg=Bg)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=as.numeric(fish$estimate),
                         down=fish$conf.int[1],
                         up=fish$conf.int[2])
      res0
   })
###
  df$odds <- res$odds
  df$down <- res$down
  df$up <- res$up  
  df  
}

ExampleGOplot <- function(cg, nbreak){

### prepare data    
   ## x <- str_split(cg$GeneRatio, "/", simplify=T)
   ## GeneRatio <- as.numeric(x[,1])/as.numeric(x[,2])
   ## Drt2 <- c("Up"=1, "Down"=2) 
   ## cg <- cg%>%mutate(Direction2=Drt2[direction],
   ##    contrast2=paste(Direction2, contrast, sep="."))%>%
   ##    mutate(contrast2=gsub("-", "+", contrast2)) 
   ## cg$size <- rep(1,nrow(cg))
   ## cg$size[GeneRatio>=0.05&GeneRatio<0.15] <- 2
   ## cg$size[GeneRatio>=0.15] <- 3 
   #
   ## cg <- cg%>%drop_na(odds)
   ## cg <- cg%>%mutate(p2=ifelse(p.adjust<pth, p.adjust, NA))
    
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls))+
      geom_point(aes(size=odds, colour=p2))+
      scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
         "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
         "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
         "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
       ## scale_colour_gradient(name="p.adjust",
       ##     low="blue", high="red", na.value=NA, trans="reverse",
       ##     n.breaks=5,           
       ##     guide=guide_colourbar(barwidth=grid::unit(0.4,"lines"),
       ##        barheight=grid::unit(1.5,"lines"), order=1))+    #"#ffa500"
      ## scale_colour_gradientn(name="p.adjust",
      ##      colours=mycolor, na.value=NA, trans="reverse",
      ##      guide=guide_colourbar(barwidth=grid::unit(0.4,"lines"),
      ##         barheight=grid::unit(2,"lines"),
      ##         order=1))+    #"#ffa500"       
       scale_size_binned("odds ratio",
                         breaks=waiver(), n.breaks=nbreak,
                         guide=guide_bins(show.limits=TRUE, axis=TRUE,
                             axis.show=arrow(length=unit(1.5,"mm"), ends="both"),
                             keywidth=grid::unit(0.4,"lines"),
                             keyheight=grid::unit(0.4,"lines"),order=2))+
        theme_bw()
    ## +
         ## theme(axis.title=element_blank(),
         ##       legend.background=element_blank(),
         ##       legend.title=element_text(size=8),
         ##       legend.text=element_text(size=7),
         ##       legend.key.size=grid::unit(0.4, "lines"))
         ##     ## legend.key.size=grid::unit(0.6, "lines"))
   fig0
}

###
###
getdata <- function(cg, tmp){
   ## cg2 <- cg%>%filter(Description==pathway)
   cg2 <- cg 
   cg2 <- cg2%>%as.data.frame()%>%
      mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
             Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
             Diff.not=Diff.total-Diff.in)
   cg2 <- cg2%>%
       mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
              Bg.total=as.numeric(gsub(".*/","", BgRatio)),
              Bg.not=Bg.total-Bg.in)

   cg2 <- odds.fun(cg2)
   cg2 <- cg2%>%dplyr::select(Cluster, p.adjust, odds)%>%full_join(tmp, by=c("Cluster"="rn"))
   ## 
   Drt2 <- c("Up"=1, "Down"=2) 
   cg2 <- cg2%>%mutate(Direction2=Drt2[direction],
      contrast2=paste(Direction2, contrast, sep="."))%>%
      mutate(contrast2=gsub("-", "+", contrast2))  
   ## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
   ## cg2$p2 <- cg2$p.adjust
   ## cg2$p2[cg2$p2>0.1] <- NA
   cg2
}

###
 
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
rn <- paste(rep(contrast, each=8), rep(rep(MCls, each=2), times=4),
            rep(rep(c("Down", "Up"),times=4), times=4), sep=".")
tmp <- data.frame(contrast=rep(contrast, each=8),
   MCls=rep(rep(MCls, each=2), times=4),
   direction=rep(rep(c("Down", "Up"),times=4), times=4))%>%
   mutate(rn=paste(contrast, MCls, direction, sep="."))

###
###
cg <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/1_enrichGO.rds")
cg <- as.data.frame(cg)
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")
ck <- as.data.frame(ck)
### shared across cell types and treatments


### Type I interferon  ###
df.path <- data.frame(pathway=c("type I interferon signaling pathway",
   "Coronavirus disease - COVID-19",
   "cellular response to glucocorticoid stimulus","Asthma"),
   labels=c("Type I IFN signaling",
            "COVID-19",
            "cellular response to glucocorticoid",
            "Asthma"))
   ## size=c(10,10,10,10),
   ## axis.x=c(F,F,T,T),
   ## breaks=c(3,2,2,2),
   ## pths=c(0.05, 0.2, 0.5, 0.5))

###
###

dots <- list()


### Type I IFN 
i <- 1
pathway0 <- df.path[i,1]
label0 <- df.path[i,2]
###
cg2 <- cg%>%filter(Description==pathway0)
cg2 <- getdata(cg2, tmp=tmp)
##
mycolor <- colorRampPalette(c("blue", "red"))(3)
pp <- cg2$p.adjust
p2 <- rep(NA, length(pp))
p2[pp<=0.2] <- 1
p2[pp<=0.1] <- 2
p2[pp<0.01] <- 3
cg2$p2 <- as.character(p2) 

p0 <- ExampleGOplot(cg2, nbreak=3)+
   scale_colour_manual(name="p.adjust",
      values=c("1"=mycolor[1], "2"=mycolor[2], "3"=mycolor[3]), 
      labels=c("1"="0.1~0.2", "2"="0.01~0.1", "3"="~0.01"), na.value=NA,
      guide=guide_legend(override.aes=list(size=2), order=1))+
   ggtitle(label0)+ 
   theme(axis.title=element_blank(),
         axis.text.y=element_text(size=9),
         axis.text.x=element_text(angle=-45, size=9, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=10),
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.4, "lines"))
dots[[i]] <- p0 



####
#### covid-19
i <- 2
pathway0 <- df.path[i,1]
label0 <- df.path[i,2]
###
cg2 <- ck%>%filter(Description==pathway0)
cg2 <- getdata(cg2, tmp=tmp)
##
mycolor <- colorRampPalette(c("blue", "red"))(3)
pp <- cg2$p.adjust
p2 <- rep(NA, length(pp))
p2[pp<=0.21] <- "1"
p2[pp<=0.1] <- "2"
p2[pp<1e-3] <- "3"
cg2$p2 <- as.character(p2) 
 
p0 <- ExampleGOplot(cg2, nbreak=3)+
   scale_colour_manual(name="p.adjust",
      values=c("1"=mycolor[1], "2"=mycolor[2], "3"=mycolor[3]), 
      labels=c("1"="0.1~0.2", "2"="1e-3~0.1", "3"="~1e-3"), na.value=NA,
      guide=guide_legend(override.aes=list(size=2), order=1))+
   ggtitle(label0)+ 
   theme(axis.title=element_blank(),
         axis.text.y=element_text(size=9),
         axis.text.x=element_text(angle=-45, size=9, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=10),
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.4, "lines"))
dots[[i]] <- p0


###
### response to glucocorticoid
i <- 3
pathway0 <- df.path[i,1]
label0 <- df.path[i,2]
###
cg2 <- cg%>%filter(Description==pathway0)
cg2 <- getdata(cg2, tmp=tmp)
##
mycolor <- colorRampPalette(c("blue", "red"))(3)
pp <- cg2$p.adjust
p2 <- rep(NA, length(pp))
p2[pp<=0.2] <- "1"
p2[pp<=0.1] <- "2"
p2[pp<0.05] <- "3"
cg2$p2 <- as.character(p2) 
 
p0 <- ExampleGOplot(cg2, nbreak=3)+
   scale_colour_manual(name="p.adjust",
      values=c("1"=mycolor[1], "2"=mycolor[2], "3"=mycolor[3]),
      labels=c("1"="0.1~0.2", "2"="0.05~0.1", "3"="~0.05"), na.value=NA,
      guide=guide_legend(override.aes=list(size=2), order=1))+
   ggtitle(label0)+ 
   theme(axis.title=element_blank(),
         axis.text.y=element_text(size=9),
         axis.text.x=element_text(angle=-45, size=9, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=10),
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.4, "lines"))
dots[[i]] <- p0


###
###
i <- 4
pathway0 <- df.path[i,1]
label0 <- df.path[i,2]
###
cg2 <- ck%>%filter(Description==pathway0)
cg2 <- getdata(cg2, tmp=tmp)
##
mycolor <- colorRampPalette(c("blue", "red"))(3)
pp <- cg2$p.adjust
p2 <- rep(NA, length(pp))
p2[pp<=0.2] <- "1"
p2[pp<=0.1] <- "2"
p2[pp<0.05] <- "3"
cg2$p2 <- as.character(p2) 

p0 <- ExampleGOplot(cg2, nbreak=3)+
   scale_colour_manual(name="p.adjust",
      values=c("1"=mycolor[1], "2"=mycolor[2], "3"=mycolor[3]),
      labels=c("1"="0.1~0.2", "2"="0.05~0.1", "3"="~0.05"), na.value=NA,
      guide=guide_legend(override.aes=list(size=2), order=1))+
   ggtitle(label0)+ 
   theme(axis.title=element_blank(),
         axis.text.y=element_text(size=9),
         axis.text.x=element_text(angle=-45, size=9, hjust=0, vjust=0.5),
         plot.title=element_text(hjust=0.5, size=10),
         legend.background=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=6),
         legend.key.size=grid::unit(0.4, "lines"))
dots[[i]] <- p0







#####################
### pathway score ###
#####################
###
avePathway <- function(X){
### filtering more missing value   
   ii <- apply(!is.na(X), 1, sum)
   X <- X[ii>20,]
   bti <- colnames(X)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1],treats=cvt[,2],sampleID=cvt[,3],Batch=cvt[,4])%>%mutate(comb=paste(MCls, Batch, sep="_"))
   comb <- unique(cvt$comb)
   for (ii in comb){
      bti0 <- cvt%>%filter(comb==ii)%>%dplyr::pull(rn)
      x <- X[,bti0]
      x.mean <- apply(x, 1, mean, na.rm=T)
      x.scale <- sweep(x, 1, x.mean, "-")
      X[,bti0] <- x.scale
      ###ctrl
      x <- X[,bti0]
      bti2 <- cvt%>%filter(comb==ii,treats=="CTRL")%>%dplyr::pull(rn)
      x0 <- apply(X[,bti2], 1, mean, na.rm=T)
      x.scale <- sweep(x, 1, x0, "-")
      X[,bti0] <- x.scale
   }
    
   pathway <- apply(X, 2, mean, na.rm=T)
}

###bulk, NB.mu, NB.phi
getPathwayData <- function(gene){

### bulk data
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
###
   cvt
}

###
###
df.path <- data.frame(pathway=c("type I interferon signaling pathway",
   "Coronavirus disease - COVID-19",
   "cellular response to glucocorticoid stimulus","Asthma"),
   labels=c("Type I IFN signaling",
            "COVID-19",
            "cellular response to glucocorticoid",
            "Asthma"),
   size=c(10,10,10,10),axis.x=c(F,F,T,T))

scores <- lapply(1:nrow(df.path), function(i){
###
    pathway0 <- df.path[i,1]
    label0 <- df.path[i,2]
    size0 <- df.path[i,3]
### 
   cg2 <- cg%>%filter(Description==pathway0)%>%as.data.frame()
   if( nrow(cg2)==0) cg2 <- ck%>%filter(Description==pathway0)%>%as.data.frame()  
   geneList <- unique(unlist(str_split(cg2$geneID,"/")))    
   gene2 <- bitr(geneList, fromType="ENTREZID", toType=c("ENSEMBL", "SYMBOL"), OrgDb=org.Hs.eg.db)
   ens <- gene2%>%dplyr::pull(ENSEMBL)
###    
   cvt <- getPathwayData(ens)%>%
       mutate(treat2=gsub("-EtOH", "", treats))    
                          
   cvt2 <- cvt%>%filter(treat2!="CTRL")#%>%drop_na(y)
   p0 <- ggplot(cvt2, aes(x=MCls, y=y, fill=treat2))+
      geom_boxplot(outlier.shape=NA)+ 
      scale_fill_manual("",
         values=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
               "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"),
         labels=c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
               "PHA"="PHA", "PHA-DEX"="PHA+DEX"))+
      ## scale_colour_manual(values=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
      ##          "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+ 
      ylab("Pathway score")+      
      ## ggtitle(label0)+
      theme_bw()+
      theme(legend.position="none",
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=10),
         axis.text.x=element_text(size=10),
         axis.text.y=element_text(size=9),
         plot.title=element_blank())
    ###
    ## if (i==1){
    ##    p0 <- p0+annotate("text", x=2.2, y=1.75,
    ##   label=bquote(atop(italic(JAK1)~","~italic(STAT1)~","~italic(STAT2)~","
    ##   ~italic(IRF1)~","~italic(IRF7)~",", italic(IRF8)~","~italic(MX1)~","
    ##   ~italic(OAS1)~" ...")), size=2.5, parse=FALSE)   
    ## }        
    p0
})



####

## figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.2_pathway.pdf"
## fig2 <- plot_grid(dots[[1]], scores[[1]],
##                   dots[[2]], scores[[2]],
##                   nrow=1, ncol=4,
##                   align="h", axis="tb",
##                   rel_widths=c(1.2, 1, 1.2, 1),
##                   labels=c("C", "D", "E", "F"),
##                   label_fontface="plain")

## fig3 <- plot_grid(dots[[3]], scores[[3]],
##                   dots[[4]], scores[[4]],
##                   nrow=1, ncol=4,
##                   align="h", axis="tb",
##                   rel_widths=c(1.2, 1, 1.2, 1),
##                   labels=c("G", "H", "I", "J"),
##                   label_fontface="plain")
## ###
## pdf(figfn, width=12, height=5)
## plot_grid(fig2, fig3, nrow=2,
##           align="v", rel_heights=c(1,1.2), labels=NULL)
## dev.off()


###
### dots[[1]]
fig1 <- plot_grid(dots[[1]], scores[[1]],
                  nrow=2, ncol=1, align="v", axis="lr",
                  rel_heights=c(1.2,1.1))
                  ## label_size=16,
                  ## label_x=c(0.1,0.1),label_y=c(1,1.1),
                  ## labels=NULL, label_fontface="plain")
fig2 <- plot_grid(dots[[2]], scores[[2]],
                  nrow=2, ncol=1, align="v", axis="lr",
                  rel_heights=c(1.2,1.1))
                  ## label_size=16,
                  ## label_x=c(0.1,0.1),label_y=c(1,1.1),
                  ## labels=NULL, label_fontface="plain")
fig3 <- plot_grid(dots[[3]], scores[[3]],
                  nrow=2, ncol=1, align="v", axis="lr",
                  rel_heights=c(1.2,1.1))
                  ## label_size=16,
                  ## label_x=c(0.08,0.08), label_y=c(1,1.1),
                  ## labels=NULL, label_fontface="plain")
fig4 <- plot_grid(dots[[4]], scores[[4]],
                  nrow=2, ncol=1, align="v", axis="lr",
                  rel_heights=c(1.2,1.1))
                  ## label_size=16,
                  ## label_x=c(0.1,0.1),label_y=c(1, 1.1),
                  ## labels=NULL, label_fontface="plain")


## figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.3_pathway.png"
## png(figfn, width=1200, height=500, res=120)
## plot_grid(fig2, fig3, fig4, fig5,
##           nrow=1, ncol=4,
##           align="h", axis="b",
##           rel_widths=1, labels=NULL)
## dev.off()

figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.3.1_pathway_review.png"
png(figfn, width=480, height=620, res=120)
print(fig1)
dev.off()

###
figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.3.2_pathway_review.png"
png(figfn, width=480, height=620, res=120)
print(fig2)
dev.off()

###
figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.3.3_pathway_review.png"
png(figfn, width=480, height=620, res=120)
print(fig3)
dev.off()

##
figfn <- "./6_DEG.CelltypeNew_output/Filter2_pub/Figure2.3.4_pathway_review.png"
png(figfn, width=480, height=620, res=120)
print(fig4)
dev.off()

