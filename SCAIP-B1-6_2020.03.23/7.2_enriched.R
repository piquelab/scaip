###
###
library(tidyverse)
library(purrr)
library(furrr) 
library("BiocParallel")
library(qqman)
library(qvalue)
##
## library(DESeq2)
## library(annotables)
## library(biobroom)
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
## theme_set(theme_grey())

rm(list=ls())

outdir <- "./7_GSE.ClusterProfiler_output/Filter2_New/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


#############################################
### supplement figures, specific GO terms ###
#############################################

###
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

ExampleGOplot <- function(cg){

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
    
   contrast2_value <- c("LPS.Up"=1, "PHA.Up"=2, "LPS+DEX.Up"=3, "PHA+DEX.Up"=4,
                        "LPS.Down"=5, "PHA.Down"=6, "LPS+DEX.Down"=7, "PHA+DEX.Down"=8) 
   cg <- cg%>%mutate(contrast2=gsub("-", "+", paste(contrast, direction, sep=".")),
                       contrast2_value=contrast2_value[contrast2],
                       contrast2=fct_reorder(contrast2, contrast2_value))
    
    
   fig0 <- ggplot(cg, aes(x=contrast2, y=MCls))+
      geom_point(aes(size=odds, colour=p2))+
      ## scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
      ##    "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
      ##    "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
      ##    "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      ## scale_colour_manual(name="p.adjust", values=c("1"="blue", "2"="red"),
      ##       labels=c("1"=">0.1", "2"="<0.1"), guide=guide_legend(override.aes=list(size=2)))+ 
       scale_colour_gradient(name="p.adjust",
           low="blue", high="red", na.value=NA, trans="reverse",
           breaks=waiver(), n.breaks=5,
           guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
              barheight=grid::unit(4,"lines"),
              order=1))+    #"#ffa500"
       scale_size_binned("odds ratio",
            breaks=waiver(), n.breaks=3,
            guide=guide_bins(show.limits=TRUE, axis=TRUE,
                axis.show=arrow(length=unit(1.5,"mm"), ends="both"),
                keywidth=grid::unit(0.6,"lines"),
                keyheight=grid::unit(0.6,"lines"),order=2))+    
      ## scale_colour_gradient(name=bquote(~log[2]~"(Odds ratio)"),                           
      ##    low="blue", high="red", na.value=NA, n.breaks=5,
      ##    guide=guide_colourbar(order=1))+    #"#ffa500"
      ## scale_size_manual(name=bquote(~log[10]~"("~P[adj]~")"),
      ##    values=c("1"=1, "2"=2.5, "3"=3.5),
      ##    labels=c("1"="<1", "2"="1~2", "3"=">3"), 
      ##    guide=guide_legend(order=2))+
      theme_bw()+
      theme(axis.title=element_blank(),
         axis.text.y=element_text(size=12),
         legend.background=element_blank(),
         legend.title=element_text(size=10),
         legend.text=element_text(size=8))
         ## legend.key.size=grid::unit(0.6, "lines"))
   fig0
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
ck <- read_rds("./7_GSE.ClusterProfiler_output/Filter2/2_enrichKEGG.rds")


###
###
getdata <- function(cg, pathway, tmp){
   cg2 <- cg%>%filter(Description==pathway)
   cg2 <- cg2%>%as.data.frame()%>%
      mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
             Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
             Diff.not=Diff.total-Diff.in)
   cg2 <- cg2%>%
       mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
              Bg.total=as.numeric(gsub(".*/","", BgRatio)),
              Bg.not=Bg.total-Bg.in)

   cg2 <- odds.fun(cg2)
   cg2 <- cg2%>%dplyr::select(Cluster, Diff.in, Diff.total, Diff.not,
       Bg.in, Bg.total, Bg.not, p.adjust, odds)%>%full_join(tmp, by=c("Cluster"="rn"))
   ###
   cg2 <- cg2%>%mutate(p2=ifelse(p.adjust<0.2, p.adjust, NA))
###  
   cg2
}



################################################
### combine specific GO terms plots together ###
################################################

pathway_ls <- c("response to lipopolysaccharide",
   "cytokine-mediated signaling pathway",
   "innate immune response",
   "stress response to metal ion")

fig_ls <- lapply(pathway_ls, function(pathway){

  if (pathway=="Asthma"){
     cg2 <- getdata(ck, pathway, tmp=tmp)
  }else{
     cg2 <- getdata(cg, pathway, tmp=tmp)
  }   
                          
  ##
  fig <- ExampleGOplot(cg2)+
     ggtitle(pathway)+      
     theme(plot.title=element_text(hjust=0.5, size=12),
           axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=10))
   fig
})

#### supplementary file
figfn <- "./7_GSE.ClusterProfiler_output/Filter2_New/Figure5_comb_response_review.pdf"
pdf(figfn, width=10, height=8, pointsize=20)
print(plot_grid(plotlist=fig_ls, nrow=2, ncol=2, labels="AUTO", label_fontface="plain"))
dev.off()



##########################################
### DEX-specific, response to reviwers ###
##########################################

pathway0 <- "stress response to metal ion"
cg2 <- getdata(cg, pathway=pathway0, tmp=tmp)
##
p <- ExampleGOplot(cg2)+
   ggtitle(pathway0)+
   theme(plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=10))

figfn <-"./7_GSE.ClusterProfiler_output/Filter2_New/Figure5_review.pdf"
pdf(figfn, width=7, height=5)
print(p)
dev.off()






### shared across cell types and treatments

##########################
### Type I interferon  ###
##########################

## cg2 <- getdata(cg, pathway="type I interferon signaling pathway", tmp=tmp)

## p1 <- ExampleGOplot(cg2)+
##    ggtitle("Type I IFN signaling")+
##    theme(axis.text.x=element_blank(),
##          plot.title=element_text(hjust=0.5, size=12))

## figfn <- paste("./7_GSE.ClusterProfiler_output/Filter2_New/Figure1.1_",
##    pathway, ".png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig)
## dev.off()

########################
### COVID-19 pathway ###
########################

## pathway <- 
## ck2 <- getdata(ck, pathway="Coronavirus disease - COVID-19", tmp=tmp)
## p2 <- ExampleGOplot(ck2)+
##    ggtitle("COVID-19")+
##    theme(axis.text.x=element_text(angle=-90, size=8, hjust=0, vjust=0.5),
##          plot.title=element_text(hjust=0.5, size=12))

## figfn <- "./7_GSE.ClusterProfiler_output/Filter2_New/Figure1.0_IFNandCOVID19.png"
## png(figfn, width=500, height=400, res=120)
## print(plot_grid(p1, p2, nrow=2, ncol=1, labels=C("C", "E"), label_fontface="plain"))
## dev.off()


#### DEX-specific 
#### stress response to metal ion
## cg2 <- cg%>%filter(Description=="stress response to metal ion")

## cg2 <- cg2%>%as.data.frame()%>%
##     mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##            Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##            Diff.not=Diff.total-Diff.in)
## cg2 <- cg2%>%
##     mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##            Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##            Bg.not=Bg.total-Bg.in)

## cg2 <- odds.fun(cg2)
## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))

## fig <- ExampleGOplot(cg2)+
##    ggtitle("stress response to metal ion")+
##    theme(plot.title=element_text(hjust=0.5, size=12))

## figfn <-"./7_GSE.ClusterProfiler_output/Filter2_pub/Figure3_metal.ion.png"
## png(figfn, width=500, height=400, res=120)
## print(fig)
## dev.off()

## ###
## ### glucocorticoid
## pathway <- "cellular response to glucocorticoid stimulus"
## cg2 <- cg%>%filter(Description==pathway)

## cg2 <- cg2%>%as.data.frame()%>%
##     mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##            Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##            Diff.not=Diff.total-Diff.in)
## cg2 <- cg2%>%
##     mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##            Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##            Bg.not=Bg.total-Bg.in)

## cg2 <- odds.fun(cg2)
## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig <- ExampleGOplot(cg2)+
##    ggtitle(pathway)+
##    theme(plot.title=element_text(hjust=0.5, size=11))

## figfn <- paste("./7_GSE.ClusterProfiler_output/Filter2_New/Figure2.2_",
##    pathway, ".png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig)
## dev.off()


## ###
## ### cell type specific
## ### leukocyte degranulation
## ### myeloid cell activation involved in immune response
## ### myeloid leukocyte mediated immunity
## ### granulocyte activation
## ### neutrophil degranulation
## ### neutrophil activation involved in immune response
## ### neutrophil activation
## pathway <- "neutrophil activation"
## i <- 7
## cg2 <- cg%>%
##    filter(Description==pathway)

## cg2 <- cg2%>%as.data.frame()%>%
##     mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##            Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##            Diff.not=Diff.total-Diff.in)
## cg2 <- cg2%>%
##     mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##            Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##            Bg.not=Bg.total-Bg.in)

## cg2 <- odds.fun(cg2)
## cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
## cg2$p2 <- cg2$p.adjust
## cg2$p2[cg2$p2>0.1] <- NA

## fig <- ExampleGOplot(cg2)+
##    ggtitle(pathway)+
##    theme(plot.title=element_text(hjust=0.5, size=12))

## figfn <-paste("./7_GSE.ClusterProfiler_output/Filter2_New/Figure3.", i,"_",
##    pathway, ".png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig)
## dev.off()

## ###
## pathway <- "Asthma"
## i <- 8
## ck2 <- ck%>%
##    filter(Description==pathway)

## ck2 <- ck2%>%as.data.frame()%>%
##     mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
##            Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
##            Diff.not=Diff.total-Diff.in)
## ck2 <- ck2%>%
##     mutate(Bg.in=as.numeric(gsub("/.*","", BgRatio)),
##            Bg.total=as.numeric(gsub(".*/","", BgRatio)),
##            Bg.not=Bg.total-Bg.in)

## ck2 <- odds.fun(ck2)
## ck2 <- ck2%>%full_join(tmp, by=c("Cluster"="rn"))
## ck2$p2 <- ck2$p.adjust
## ck2$p2[ck2$p2>0.1] <- NA

## fig <- ExampleGOplot(ck2)+
##    ggtitle(pathway)+
##    theme(plot.title=element_text(hjust=0.5, size=12))

## figfn <-paste("./7_GSE.ClusterProfiler_output/Filter2_New/Figure3.", i,"_",
##    pathway, ".png", sep="")
## png(figfn, width=500, height=400, res=120)
## print(fig)
## dev.off()
