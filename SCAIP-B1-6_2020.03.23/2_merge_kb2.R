#
source("./Bin/LibraryPackage.R")

rm(list=ls())
outdir <- "./2_kb2_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

## read data including barcodes.txt, genes.txt and mtx
## 2_merge_kb2 for downstream analysis

################################################
### 1, generate folders containing h5ad data ###
################################################
if(FALSE){
cat("1.", "generate folders", "\n")
basefolder <- "/nfs/rprdata/julong/SCAIP/kallisto2/bus/"
#basefolder <- "/nfs/rprdata/scaip/kallisto/bus/"
#basefolder <- "/nfs/rprdata/scaip/kallisto2/bus/"
expNames <- dir(basefolder,"^SCAIP*")
folders <- paste0(basefolder, expNames, "", sep="")
ind <- dir.exists(folders) #ind <- file.info(folders)$isdir;ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
#cat(folders, sep="\n")
## remove bad emulsions
EXP_del <- c("SCAIP2-LPS-DEX", "SCAIP2-PHA-DEX", "SCAIP2-PHA-EtOH", "SCAIP3-LPS-DEX", "SCAIP3-PHA-DEX", "SCAIP3-PHA-EtOH")
folders <- folders[!expNames%in%EXP_del]
expNames <- names(folders)
} ###1, End


###########################################################
### 2, read h5ad data into seurat then merge 39 objects ###
###########################################################

#library(furrr)
#future::plan(strategy = 'multicore', workers = 5)
#options(future.globals.maxSize = 10 * 1024 ^ 3)

if(FALSE){
cat("2.1.", "Read data", "\n")
## 2.2, read mtx data into seurat object
## Function to read kallisto?
readKallisto  <- function (run, prefixFile="spliced/s",expPrefix=NULL) 
{
    if (!dir.exists(paths = run)) {
        stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, paste0(prefixFile,".barcodes.txt"))
    gene.loc <- file.path(run, paste0(prefixFile,".genes.txt"))
##    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, paste0(prefixFile,".mtx"))
    if (!file.exists(barcode.loc)) {
        stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)) {
        stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
        stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if(is.null(expPrefix)){
        rownames(data) <- cell.names
    }
    else{
        rownames(data) <- paste0(expPrefix,"_",cell.names)
    }    
    feature.names <- readLines(gene.loc)
    colnames(data) <- feature.names
    t(data)
}

###
adata <- future_map(expNames,function(ii){
    ##
    expPrefix = ii;
    cat("#Loading ",paste0(folders[ii], "/counts_unfiltered/spliced"), " ...")    
    sFull <- readKallisto(folders[ii], prefixFile="/counts_unfiltered/spliced" , expPrefix)  #"/spliced/s"
    cat(dim(sFull),"\n")
    cat("#Loading ",paste0(folders[ii], "/counts_unfiltered/unspliced"), " ...")
    uFull <- readKallisto(folders[ii], prefixFile="/counts_unfiltered/unspliced", expPrefix) #"/unspliced/u" 
    cat(dim(uFull),"\n")
    ##
    scs <- colSums(sFull)
    rownames(sFull) <- paste0("S-",rownames(sFull))
    ucs <- colSums(uFull)
    rownames(uFull) <- paste0("U-",rownames(uFull))
    
    sel <- intersect(colnames(sFull),colnames(uFull))
    #sel <- intersect(colnames(sFull)[scs>0],colnames(uFull)[ucs>0])
    
    count0 <- rbind(sFull[,sel],uFull[,sel]) 
    #sc0 <- CreateSeuratObject(count0)
    #sc0
    ## May need to rename rownames or split...
    sce <- create_SCE(count0)
    cat(dim(sce),"\n")
    ## Remove debris...
    sce <- diem(sce,top_n = 16000)
    sc <- convert_to_seurat(sce)
    cat("#Final: ",dim(sc),"\n")
    sc
})

sc <- merge(adata[[1]],adata[-1], project="kbSCAIP2")
##  
###       
library(annotables)
anno <- tibble(rn=rownames(sc)) %>% mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), uns=grepl("U-",rn)) %>% left_join(grch38)                  
sc[["percent.mt"]] <- PercentageFeatureSet(sc, features = anno %>% filter(chr=="MT") %>% dplyr::pull(rn) )
sc@meta.data$NEW_BARCODE <- colnames(sc)
            
write_rds(sc, "./2_kb2_output/1_Seurat_kb.rds")

} ###End, 2.1

### 1_Seurat_kb.rds, unfiltered data and removing diem ## default data, 301,637 barcodes
### 1_Seurat_kb2.rds, unfilered data and removing dime, scs>0 and ucs>0, 304,360 barcodes
### 1_Seurat_kb3.rds, filtered data ## filtered by bustools, 306,299 barcodes 


###########################################
### 2.2, show summary stats of raw data ###
###########################################
###(1)
if (TRUE){
cat("2.2.", "Summary of raw data", "\n")

sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
count <- sc@assays$RNA@counts
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("S-",rn),rnz=rowSums(count))
geneSel <- anno%>%filter(uns)%>%dplyr::select(rn)%>%unlist()           
Xs <- count[geneSel,]
nCount_spliced <- colSums(Xs)
XsN <- Xs>0
nFeature_spliced <- colSums(XsN)           

##number of genes
tmp <- anno%>%filter(uns,rnz>0)

meta <- sc@meta.data
meta$nCount_spliced <- nCount_spliced
meta$nFeature_spliced <- nFeature_spliced

dd <- meta%>%group_by(orig.ident)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       S_reads=mean(nCount_spliced),
                       S_ngene=mean(nFeature_spliced), .groups="drop")
dd <- dd%>%dplyr::rename(ident=orig.ident)%>%
           mutate(batch=gsub("-.*","",ident))

              
###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("#Barcodes per experiment")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png("./2_kb2_output/Figure1.1_barcodes.png", width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="#UMIs per cell", "2"="#Genes per cell"))
fig0 <- ggplot(ddnew, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
              strip.background=element_blank())
##        
png("./2_kb2_output/Figure1.2_genes.png", width=1200, height=1000, res=150)
print(fig0)
dev.off() 


###(3), spliced reads and genes          
dd1 <- dd%>%dplyr::select(ident, S_reads, batch)%>%mutate(stats=1)%>%dplyr::rename(y=S_reads)
dd2 <- dd%>%dplyr::select(ident, S_ngene, batch)%>%mutate(stats=2)%>%dplyr::rename(y=S_ngene)
ddnew2 <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="#UMIs per cell", "2"="#Genes per cell"))
fig0 <- ggplot(ddnew2, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+  
        #guides(fill=guide_legend(override.aes=list(size=2)))+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=7),
              strip.background=element_blank())
###
###
png("./2_kb2_output/Figure1.3_spliced.png", width=1200, height=1000, res=150)
print(fig0)
dev.off()

     
### (4), unspliced genes  ###
rm(list=ls())
sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
count <- sc@assays$RNA@counts
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("U-",rn),
               rnz=rowSums(count))
geneSel <- anno%>%filter(uns,rnz>0)%>%dplyr::select(rn)%>%unlist()           
Xu <- count[geneSel,]
nCount_unspliced <- colSums(Xu)
nFeature_unspliced <- colSums(Xu>0)           

tmp <- anno%>%filter(uns,rnz>0) ##24,526 genes

meta <- sc@meta.data
meta$nCount_unspliced <- nCount_unspliced
meta$nFeature_unspliced <- nFeature_unspliced

dd2 <- meta%>%group_by(orig.ident)%>%
             summarise(ncell=n(),
                       U_reads=mean(nCount_unspliced),
                       U_ngene=mean(nFeature_unspliced),.groups="drop")
dd2 <- dd2%>%dplyr::rename(ident=orig.ident)%>%
           mutate(batch=gsub("-.*","",ident))

dd2a <- dd2%>%dplyr::select(ident, U_reads, batch)%>%mutate(stats=1)%>%dplyr::rename(y=U_reads)
dd2b <- dd2%>%dplyr::select(ident, U_ngene, batch)%>%mutate(stats=2)%>%dplyr::rename(y=U_ngene)
ddnew2 <- rbind(dd2a,dd2b)      

stats <- as_labeller(c("1"="#UMIs per cell","2"="#Genes per cell"))
fig0 <- ggplot(ddnew2, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+  
        #guides(fill=guide_legend(override.aes=list(size=2)))+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=7),
              strip.background=element_blank())
###
###
png("./2_kb2_output/Figure1.4_unspliced.png", width=1200, height=1000, res=150)
print(fig0)
dev.off() 
} ## End, 2.2

     

#######################################
### 3, filter data by demux results ###
#######################################
if(FALSE){
rm(list=ls())

##########################
### 3.1, filtered data ###
##########################
cat("3.1.", "filter data by removing mismatching barcodes", "\n")

sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
meta <- sc@meta.data
meta$NEW_BARCODE <- rownames(meta)
demux <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
meta <- meta%>%inner_join(demux,by="NEW_BARCODE")

x <- read.table("../CorrectBatch1-6.txt", header=T)
correctBatch <- paste(x$sample, x$Group, sep="_")

meta0 <- meta%>%
         mutate(comb=paste(BEST.GUESS,BATCH, sep="_"))%>%
         filter(comb%in%correctBatch,grepl("SNG",DROPLET.TYPE))
sc2 <- subset(sc, cells=meta0$NEW_BARCODE)
rownames(meta0) <- meta0$NEW_BARCODE
sc2@meta.data <- meta0
 
### 292,394 barcodes and 116734(spliced and unsplcied genes) 
opfn <- "./2_kb2_output/2_Seurat_kb.rds"  
write_rds(sc2,opfn)

sparse.size <- object.size(sc2)
print(sparse.size,units="GB")    ###12.2 GB
} ### End, 3.1


#################
### 3.2, show ###
#################
if(FALSE){
cat("3.2.", "Summary of clean data", "\n")
rm(list=ls())
sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
count <- sc@assays$RNA@counts
rn <- rownames(count)
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn), 
               uns=grepl("S-",rn),
               rnz=rowSums(count))
geneSel <- anno%>%filter(uns,rnz>0)%>%dplyr::pull(rn)           
Xs <- count[geneSel,]
nCount_spliced <- colSums(Xs)
nFeature_spliced <- colSums(Xs>0)           

tmp <- anno%>%filter(uns,rnz>0) ##46,384 genes

meta <- sc@meta.data
meta$nCount_spliced <- nCount_spliced
meta$nFeature_spliced <- nFeature_spliced

dd <- meta%>%group_by(orig.ident)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       S_reads=mean(nCount_spliced),
                       S_ngene=mean(nFeature_spliced),.groups="drop")
dd <- dd%>%dplyr::rename(ident=orig.ident)%>%
           mutate(batch=gsub("-.*","",ident))


###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("#Barcodes per experiment")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png("./2_kb2_output/Figure2.1_barcodes.png", width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="#UMIs per cell", "2"="#Genes per cell"))
fig0 <- ggplot(ddnew, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
              strip.background=element_blank())
###
###
png("./2_kb2_output/Figure2.2_genes.png", width=1200, height=1000, res=150)
print(fig0)
dev.off() 

###(3), spliced reads and genes          
dd1 <- dd%>%dplyr::select(ident, S_reads, batch)%>%mutate(stats=1)%>%dplyr::rename(y=S_reads)
dd2 <- dd%>%dplyr::select(ident, S_ngene, batch)%>%mutate(stats=2)%>%dplyr::rename(y=S_ngene)
ddnew2 <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="#UMIs per cell", "2"="#Genes per cell"))
fig0 <- ggplot(ddnew2, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+  
        #guides(fill=guide_legend(override.aes=list(size=2)))+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=7),
              strip.background=element_blank())
###
###
png("./2_kb2_output/Figure2.3_spliced.png", width=1200, height=1000, res=150)
print(fig0)
dev.off()

     
### (4), unspliced genes  ###
rm(list=ls())
sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
count <- sc@assays$RNA@counts
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("U-",rn),
               rnz=rowSums(count))
geneSel <- anno%>%filter(uns,rnz>0)%>%dplyr::select(rn)%>%unlist()           
Xu <- count[geneSel,]
nCount_unspliced <- colSums(Xu)
nFeature_unspliced <- colSums(Xu>0)           

tmp <- anno%>%filter(uns,rnz>0) ##24,526 genes

meta <- sc@meta.data
meta$nCount_unspliced <- nCount_unspliced
meta$nFeature_unspliced <- nFeature_unspliced

dd2 <- meta%>%group_by(orig.ident)%>%
             summarise(ncell=n(),
                       U_reads=mean(nCount_unspliced),
                       U_ngene=mean(nFeature_unspliced))
dd2 <- dd2%>%dplyr::rename(ident=orig.ident)%>%
           mutate(batch=gsub("-.*","",ident))

dd2a <- dd2%>%dplyr::select(ident, U_reads, batch)%>%mutate(stats=1)%>%dplyr::rename(y=U_reads)
dd2b <- dd2%>%dplyr::select(ident, U_ngene, batch)%>%mutate(stats=2)%>%dplyr::rename(y=U_ngene)
ddnew2 <- rbind(dd2a,dd2b)      

stats <- as_labeller(c("1"="#UMIs per cell","2"="#Genes per cell"))
fig0 <- ggplot(ddnew2, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+  
        #guides(fill=guide_legend(override.aes=list(size=2)))+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=7),
              strip.background=element_blank())
###
###
png("./2_kb2_output/Figure2.4_unspliced.png", width=1200, height=1000, res=150)
print(fig0)
dev.off() 


### (5), Table showing summary data by individuals
meta <- sc@meta.data
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*", "", rn), 
               uns=grepl("S-",rn),
               rnz=rowSums(count))
geneSel <- anno%>%filter(uns,rnz>0)%>%dplyr::select(rn)%>%unlist()           
Xs <- count[geneSel,]
nCount_spliced <- colSums(Xs)
nFeature_spliced <- colSums(Xs>0)
meta$nCount_spliced <- nCount_spliced
meta$nFeature_spliced <- nFeature_spliced
dd2 <- meta%>%
       group_by(BEST.GUESS, treats)%>%
       summarise(ncell=n(),
                 S_reads=mean(nCount_spliced),
                 S_ngene=mean(nFeature_spliced))
dd2$treats <- gsub("-EtOH", "", dd2$treats)

tmp <- dd2%>%group_by(treats)%>%summarise(nind=n(),ncell=median(ncell),reads=median(S_reads),ngene=median(S_ngene))


col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
lab1 <- c("CTRL"="CTRL",
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")

fig1 <- ggplot(dd2,aes(x=treats, y=ncell, fill=treats))+
        geom_violin()+xlab("")+ylab("")+
        ggtitle("#Cells per individual")+
        scale_fill_manual(values=col1)+
        scale_x_discrete(labels=lab1)+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
fig2 <- ggplot(dd2,aes(x=treats, y=S_reads, fill=treats))+
        geom_violin()+xlab("")+ylab("")+
        ggtitle("#UMIs per cell")+
        scale_fill_manual(values=col1)+
        scale_x_discrete(labels=lab1)+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
              
fig3 <- ggplot(dd2,aes(x=treats, y=S_ngene, fill=treats))+
        geom_violin()+xlab("")+ylab("")+
        ggtitle("#Genes per cell")+
        scale_fill_manual(values=col1)+
        scale_x_discrete(labels=lab1)+
        theme_bw()+
        theme(legend.position="none",
              plot.title=element_text(hjust=0.5),
              axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
             
png("./2_kb2_output/Figure2.5_violin.png",width=800, height=500, res=120)
print(plot_grid(fig1, fig2, fig3, ncol=3))
dev.off()
} ##3.2, End

###
###
if (FALSE){
sc <- read_rds("./2_kb2_output/2_Seurat_kb.rds")
meta <- sc@meta.data
write.csv(meta, file="./2_kb2_output/3_kb2.meta",row.names=T)
}



############################################################################################
### 4, summary mismatching barcodes and looking for some reason that lost lots barcodes  ###
############################################################################################
#rm(list=ls())

#sc <- read_rds("./2_kb2_output/1_Seurat_kb.rds")
#dd <- read_rds("./1_demux2_output/1_demux_filter0.ALL.rds")
#cellSel <- intersect(colnames(sc),dd$NEW_BARCODE)
#sc2 <- subset(sc, cells=cellSel) ## common barcode shared between kb and demuxlet
#meta <- sc2@meta.data%>%left_join(dd, by="NEW_BARCODE")


#x <- read.table("../CorrectBatch1-6.txt", header=T)
#correctBatch <- paste(x$sample, x$Group, sep="_")

#x1 <- x%>%mutate(ii=gsub("SCAIP", "", Group))%>%
#      mutate(sample2=paste0("B", ii, "_", sample))  
#sample2 <- setNames(x1$sample2, x1$sample)
#batch2 <- setNames(x1$Group,x1$sample)


### (1) figure0.1, heatmap, all data, including matching and mismatching
#mydf <- data.frame(Sample=meta$BEST.GUESS, BATCH=meta$BATCH)%>%
#        group_by(Sample, BATCH)%>%
#        summarize(Freq=n())       
#mydf <- mydf%>%
#        group_by(Sample)%>%
#        mutate(Perc=Freq/sum(Freq)*100, 
#               Sample2=sample2[as.character(Sample)])

#fig1 <- ggplot(mydf,aes(x=as.character(Sample2), y=BATCH, fill=Perc))+
#            geom_tile()+
#            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
#            xlab("")+ylab("")+
#            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
            
#png("./2_kb_output/Figure0.1.Indi.png", width=1200, height=500, res=120)
#fig1
#dev.off()

### (2) figure0.2, barplot, mismatching barcodes
#comb <- paste(meta$BEST.GUESS, meta$BATCH, sep="_")
#meta0 <- meta[!comb%in%correctBatch,]

#dd <- meta0%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH")
#mydf0 <- dd%>%
#         group_by(orig.ident)%>%
#         summarize(ncell=n())%>%
#         rename(ident=orig.ident)%>%
#         mutate(batch=gsub("-.*","",ident))
         
#fig2 <- ggplot(mydf0,aes(x=ident, y=ncell, fill=factor(batch)))+
#        geom_bar(stat="identity")+
#        xlab("")+ylab("Number Barcodes")+
#        geom_text(aes(label=ncell),vjust=-0.5, size=3)+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              axis.text.x=element_text(angle=90, hjust=1, size=8))  
#png("./2_kb_output/Figure0.2.miss.png", width=1000, height=600, res=120)
#fig2
#dev.off()

### (3), heatmap, mismatching barcodes
#batch2 <- setNames(x1$Group,x1$sample)
#dd <- meta0%>%select("BEST.GUESS", "EXP","BATCH")%>%mutate(BATCH2=batch2[BEST.GUESS])
#mydf0 <- dd%>%group_by(EXP,BATCH2)%>%
#         summarize(ncell=n())%>%
#         group_by(EXP)%>%
#         mutate(Perc=ncell/sum(ncell)*100)

#fig3 <- ggplot(mydf0,aes(x=as.character(EXP), y=BATCH2, fill=Perc))+
#            geom_tile()+
#            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
#            xlab("")+ylab("")+
#            theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,size=6))
            
#png("./2_kb_output/Figure0.3.png", width=1000, height=500, res=120)
#fig3
#dev.off()

### (4), heatmap, mismatching barcodes
#dd <- meta0%>%select("BEST.GUESS", "EXP","BATCH")%>%mutate(sample2=sample2[BEST.GUESS])
#mydf0 <- dd%>%group_by(sample2,BATCH)%>%
#         summarize(ncell=n())%>%
#         group_by(sample2)%>%
#         mutate(Perc=ncell/sum(ncell)*100)

#fig4 <- ggplot(mydf0,aes(x=as.character(sample2), y=BATCH, fill=Perc))+
#            geom_tile()+
#            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
#            xlab("")+ylab("")+
#            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
            
#png("./2_kb_output/Figure0.4.png", width=1200, height=500, res=120)
#fig4
#dev.off()

### (5), barplot, mismatching barcodes
#comb <- paste(meta$BEST.GUESS, meta$BATCH, sep="_")
#meta0 <- meta[!comb%in%correctBatch,]

#dd <- meta0%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH")
#mydf0 <- dd%>%
#         group_by(BEST.GUESS)%>%
#         summarize(ncell=n())%>%
#         mutate(sample2=sample2[BEST.GUESS],BATCH=batch2[BEST.GUESS])
         
#fig5 <- ggplot(mydf0,aes(x=sample2, y=ncell, fill=factor(BATCH)))+
#        geom_bar(stat="identity")+
#        xlab("")+ylab("Number Barcodes")+
#        geom_text(aes(label=ncell),vjust=-0.5, size=3)+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              axis.text.x=element_text(angle=60, hjust=1, size=6))  
#png("./2_kb_output/Figure0.5.mismatching.png", width=1200, height=500, res=120)
#fig5
#dev.off()

   
### (6), barplots for matching barcodes
#comb <- paste(meta$BEST.GUESS, meta$BATCH, sep="_")
#meta1 <- meta[comb%in%correctBatch,]

#dd <- meta1%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH")
#mydf0 <- dd%>%
#         group_by(BEST.GUESS, BATCH)%>%
#         summarize(ncell=n())%>%
#         mutate(sample2=sample2[BEST.GUESS])
         
#fig6 <- ggplot(mydf0,aes(x=sample2, y=ncell, fill=factor(BATCH)))+
#        geom_bar(stat="identity")+
#        xlab("")+ylab("Number Barcodes")+
#        geom_text(aes(label=ncell),vjust=-0.5, size=3)+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              axis.text.x=element_text(angle=60, hjust=1, size=6))  
#png("./2_kb_output/Figure0.6.matching.png", width=1200, height=500, res=120)
#fig6
#dev.off()

  

### (7), barplot for all the data
#dd <- meta%>%select("orig.ident", "BEST.GUESS", "EXP", "BATCH")
#mydf0 <- dd%>%
#         group_by(BEST.GUESS)%>%
#         summarize(ncell=n())%>%
#         mutate(sample2=sample2[BEST.GUESS],BATCH=batch2[BEST.GUESS])
#         
#fig7 <- ggplot(mydf0,aes(x=sample2, y=ncell, fill=factor(BATCH)))+
#        geom_bar(stat="identity")+
#        xlab("")+ylab("Number Barcodes")+
#        geom_text(aes(label=ncell),vjust=-0.5, size=3)+
#        theme_bw()+
#        theme(legend.title=element_blank(),
#              axis.text.x=element_text(angle=60, hjust=1, size=6))  
#png("./2_kb_output/Figure0.7.all.png", width=1200, height=500, res=120)
#fig7
#dev.off()






