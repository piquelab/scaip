##
library(Matrix)
library(tidyverse)
library(data.table)
##

outdir <- "./4_summary.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=FALSE)


#################
### Read data ###
#################

## monocle
res1 <- read.csv("./1.2_monocle.outs/1_corr_pseudotime.csv")
corr1 <- res1[,2:51]
var1 <- apply(corr1, 1, var, na.rm=T)
m1 <- apply(corr1, 1, mean, na.rm=T)

## scanpy
res2 <- read.csv("./2_scanpy.outs/2_corr_pseudotime.csv")
corr2 <- res2[,2:51]
var2 <- apply(corr2, 1, var, na.rm=T)
m2 <- apply(corr2, 1, mean, na.rm=T)


## DLDA
res3 <- read.csv("./3_LDA.outs/1_corr_pseudotime.csv")
corr3 <- res3[,2:51]
var3 <- apply(corr3, 1, var, na.rm=T)
m3 <- apply(corr3, 1, mean, na.rm=T)


##########################################################
### Violin plot of variance of correlation coefficient ###
##########################################################

Df <- data.frame(gene=res1$gene, var1=var1, var2=var2, var3=var3)
nna <- rowSums(is.na(Df))
Df <- Df[nna==0,]

##
Df2 <- Df%>%pivot_longer(!gene, names_to="method", values_to="var")

###
Df2 <- Df2%>%mutate(logVar=log10(var))

##
p <- ggplot(Df2, aes(x=factor(method), y=logVar))+
    geom_violin(aes(fill=factor(method)))+
    scale_x_discrete(labels=c("var1"="Monocle", "var2"="scanpy", "var3"="DLDA"))+
    ylab(bquote(~log[10]~"(Variance of correlation)"))+
    scale_fill_manual(values=c("var1"="#7570b3", "var2"="#1b9e77", "var3"="#e7298a"))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=10))

figfn <- "./4_summary.outs/Figure1_violin.png"
png(figfn, width=380, height=420, res=120)
p
dev.off()


############################
### Violin plots of mean ###
############################

Df <- data.frame(gene=res1$gene, m1=m1, m2=m2, m3=m3)
nna <- rowSums(is.na(Df))
Df <- Df[nna==0,]

##
Df2 <- Df%>%pivot_longer(!gene, names_to="method", values_to="mm")

###
## Df2 <- Df2%>%mutate(logVar=log(var))

##
p <- ggplot(Df2, aes(x=factor(method), y=abs(mm)))+
    geom_violin(aes(fill=factor(method)))+
    scale_x_discrete(labels=c("m1"="Monocle", "m2"="scanpy", "m3"="DLDA"))+
    ylab(bquote(~"|"~rho~"|"))+
    scale_fill_manual(values=c("m1"="#7570b3", "m2"="#1b9e77", "m3"="#e7298a"))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank())

figfn <- "./4_summary.outs/Figure1.2_mean_violin.png"
png(figfn, width=420, height=420, res=120)
p
dev.off()



###############
### heatmap ###
###############

###
### top 20 genes
Df1 <- res1%>%pivot_longer(!gene, names_to="ID", values_to="rr")
x <- Df1%>%group_by(ID)%>%slice_max(abs(rr), n=20)
gene1 <- x%>%pull(gene)%>%unique()
Df1 <- Df1%>%filter(gene%in%gene1)%>%
   mutate(xvalue=as.numeric(gsub("id", "", ID)), ID=fct_reorder(ID, xvalue))

###
p <- ggplot(Df1, aes(x=factor(ID), y=factor(gene), fill=rr))+
    geom_tile()+
    xlab("Simulations")+ylab("Pseudotime-determined gene")+
    ggtitle("Monocle")+
    scale_fill_gradient2(name=bquote(rho),
       low="blue", mid="white", high="red",
       guide=guide_colorbar(barwidth=0.5, barheight=5))+
    ## coord_fixed()+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12))

###
figfn <- "./4_summary.outs/Figure2.1_monocle.heatmap.png"
png(figfn, width=320, height=420, res=120)
p
dev.off()




##scanpy
Df2 <- res2%>%pivot_longer(!gene, names_to="ID", values_to="rr")
x <- Df2%>%group_by(ID)%>%slice_max(abs(rr), n=20)
gene2 <- x%>%pull(gene)%>%unique()
Df2 <- Df2%>%filter(gene%in%gene2)%>%
   mutate(xvalue=as.numeric(gsub("id", "", ID)), ID=fct_reorder(ID, xvalue))

###
p2 <- ggplot(Df2, aes(x=factor(ID), y=factor(gene), fill=rr))+
    geom_tile()+
    xlab("Simulations")+ylab("Pseudotime-determined genes")+
    scale_fill_gradient2(name=bquote(rho),
       low="blue", mid="white", high="red",
        guide=guide_colorbar(barwidth=0.5, barheight=5))+
    ggtitle("scanpy")+
    ## coord_fixed()+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12))

###
figfn <- "./4_summary.outs/Figure2.2_scanpy.heatmap.png"
png(figfn, width=320, height=420, res=120)
p2
dev.off()



### 
Df3 <- res3%>%pivot_longer(!gene, names_to="ID", values_to="rr")
x <- Df3%>%group_by(ID)%>%slice_max(abs(rr), n=20)
gene3 <- x%>%pull(gene)%>%unique()
Df3 <- Df3%>%filter(gene%in%gene3)%>%
   mutate(xvalue=as.numeric(gsub("id", "", ID)), ID=fct_reorder(ID, xvalue))

###
p3 <- ggplot(Df3, aes(x=factor(ID), y=factor(gene), fill=rr))+
    geom_tile()+
    xlab("Simulations")+ylab("Pseudotime-determined genes")+
    scale_fill_gradient2(name=bquote(rho),
       low="blue", mid="white", high="red",
       guide=guide_colorbar(barwidth=0.5, barheight=5))+
    ggtitle("DLDA")+
    ## coord_fixed()+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=12))
 
###
figfn <- "./4_summary.outs/Figure2.3_DLDA.heatmap.png"
png(figfn, width=480, height=320, res=120)
p3
dev.off()



#####################
### Jaccard index ###
#####################

allsets <- combn(50,2)

##
cal.Jindex <- function(res, allsets, ntop){
###    
Jindex <- sapply(1:ncol(allsets), function(i){
   ##
   i1 <- allsets[1,i]+1
   r1 <- data.frame(rr=res[, i1], gene=res[,1])%>%
       drop_na(rr)%>%arrange(desc(abs(rr)))
   top1 <- r1$gene[1:ntop]
   ##
   i2 <- allsets[2,i]+1  
   r2 <- data.frame(rr=res[,i2], gene=res[,1])%>%
       drop_na(rr)%>%arrange(desc(abs(rr)))
   top2 <- r2$gene[1:ntop]
   ##
   Jindex <- length(intersect(top1, top2))/length(union(top1,top2))
   Jindex
})
###
Jindex
}


###
df1 <- data.frame(Jindex=cal.Jindex(res1, allsets, ntop=20), method="m1")
#
df2 <- data.frame(Jindex=cal.Jindex(res2, allsets, ntop=20), method="m2")
##
df3 <- data.frame(Jindex=cal.Jindex(res3, allsets, ntop=20), method="m3")

###
plotdf <- rbind(df1, df2, df3)
p <- ggplot(plotdf, aes(x=factor(method), y=Jindex))+
    geom_violin(aes(fill=factor(method)))+
    stat_summary(fun="median", geom="point", colour="red")+
    scale_x_discrete(labels=c("m1"="Monocle", "m2"="scanpy", "m3"="DLDA"))+
    ylab("Jaccard index")+
    scale_fill_manual(values=c("m1"="#7570b3", "m2"="#1b9e77", "m3"="#e7298a"))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=10))

figfn <- "./4_summary.outs/Figure3_jaccard_violin.png"
png(figfn, width=420, height=320, res=120)
p
dev.off()
    
