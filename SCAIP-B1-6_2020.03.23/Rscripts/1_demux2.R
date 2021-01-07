#
source("./Bin/LibraryPackage.R")
rm(list=ls())
outdir <- "./1_demux2_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

##
## demux2 for downstream analysis

#################################
### 1. folders of counts data ###
#################################
basefolder <- "/nfs/rprdata/julong/SCAIP/count/demuxNew/"
demuxfn <- list.files(basefolder,pattern="*.out.best")
expNames <- gsub(".out.best", "", demuxfn) 
## remove bad emulsions
#EXP_del <- c("SCAIP2-LPS-DEX", "SCAIP2-PHA-DEX", "SCAIP2-PHA-EtOH", "SCAIP3-LPS-DEX", "SCAIP3-PHA-DEX", "SCAIP3-PHA-EtOH")
#expNames <- setdiff(expNames, EXP_del)

 
############################
### 2, read demulet data ###
############################
demux <- mclapply(expNames,function(ii){
  cat("#Loading ", ii, "\n")
  fn <- paste0(basefolder, ii,".out.best")
  dd <- data.frame(fread(fn,header=T))
  dd <- dd%>%mutate(NEW_BARCODE=paste0(ii,"_", gsub("-1","",BARCODE)),EXP=ii)
},mc.cores=10)
demux <- do.call(rbind,demux)
###
demux <- demux%>%
         mutate(BEST.GUESS=gsub(",.*", "", BEST.GUESS),
                 NEXT.GUESS=gsub(",.*", "", NEXT.GUESS),
                 BATCH=substring(EXP,1,6),
                 BATCH2=gsub("-.*", "", EXP),
                 treats=gsub(".*[0-9].{,2}-","",EXP),
                 chemi=grepl("SCAIP5V3|SCAIP6",BATCH2),
                 chem=ifelse(chemi,"V3", "V2")) 
### output
opfn <- "./1_demux2_output/1_demux_New.ALL.rds" ## 1,393,290
write_rds(demux, opfn)

## 0,filtered0, AF>0&AF<1
### 1, filtered1, AF>=0.01&AF<=0.99
### 2, filtered2, AF>=0.02&AF<=0.98
### 4, not filtered 
### 5, filtered AF>0&AF<1, DP>4
### New, run samtools mpileup for 39 EXP separately, MAF>0.05&DP>100

##################################################
### 3, filter cells with incorrect information ###
##################################################

### correct batch information
#fn1 <- "/nfs/rprdata/julong/SCAIP/analyses/CorrectBatch1-5.txt"
#infor <- read.table(fn1,header=T)
#fn2 <- "/wsu/home/groups/piquelab/SCAIP/covariates/SCAIP6_dbgaps.txt"
#x <- read.table(fn2)
#batch6 <- data.frame(sample=x[,1], Group="SCAIP6")
#inforNew <- rbind(infor,batch6)

###
#opfn <- "../CorrectBatch1-6.txt"
#write.table(inforNew,file=opfn, quote=F, sep="\t", row.names=F)

##
#dd <- read_rds("./1_demux2_output/1_demux_ALL.rds") 
#ddx <- dd%>%filter(DIFF.LLK.BEST.NEXT>0)   ##267,243   
### plots
#llk <- ddx$DIFF.LLK.BEST.NEXT
#dx0 <- data.frame(llk=llk)
##
#fig0 <- ggplot(dx0,aes(x=llk))+
#        geom_density()+ggtitle("Diff.likelihood")+xlim(0,100)+
#        theme_bw()+theme(plot.title=element_text(hjust=0.5))
## 
#figfn <- "./1_demuxlet_output/figure0_DIFF.LLK.png"
#png(figfn, width=500, height=400, res=120)
#fig0
#dev.off()

dd <- read_rds("./1_demux2_output/1_demux_filter0.ALL.rds")%>%filter(DIFF.LLK.BEST.NEXT>0)
sample_batch <- paste(dd$BEST.GUESS, dd$BATCH, sep="_")
x <- read.table("../CorrectBatch1-6.txt", header=T)
correctBatch <- paste(x$sample, x$Group, sep="_")
dd0 <- dd[sample_batch %in% correctBatch,]       ###200,789

opfn <- "./1_demux2_output/2_demux_filter0.Sel.rds"
write_rds(dd0,opfn)

#opfn <- "./1_demuxlet_output/2_demuxlet_Sel2.rds"
#write_rds(ddx,opfn)



###
fn <- "./1_demux2_output/2_demux_filter0.Sel.rds"
dd <- read_rds(fn)

ddx <- dd%>%group_by(BEST.GUESS,BATCH)%>%summarise(ncell=n())%>%arrange(BATCH)
dd0 <- ddx[-17,]

fig0 <- ggplot(dd0,aes(x=BEST.GUESS, y=ncell, fill=factor(BATCH)))+
        geom_bar(stat="identity")+
        xlab("")+ylab("Number Barcodes")+
        geom_text(aes(label=ncell),vjust=-0.5, size=3)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=8))  
png("./1_demux2_output/Figure0_barcodes.png", width=1000, height=600, res=120)
fig0
dev.off()



#x <- sample_batch[sample_batch%in%correctBatch]
#myx <- str_split(x, "_", simplify=TRUE)
x <- read.table("../CorrectBatch1-6.txt", header=T)
Batch <- c("SCAIP1", "SCAIP2", "SCAIP3", "SCAIP4", "SCAIP5", "SCAIP6")
prefix <- LETTERS[1:6]
xNew <- NULL
for (i in 1:6){
xi <- x%>%filter(Group==Batch[i])
newLabel <- paste(prefix[i], xi[,1], sep="_")
tmp <- data.frame(Sample=xi[,1], Group=newLabel)
xNew <- rbind(xNew,tmp)
}
xSample2 <- setNames(xNew[,2], xNew[,1])


##mis match
mis <- sample_batch[!(sample_batch%in%correctBatch)]
mis <- str_split(mis, "_", simplify=TRUE)
mis <- data.frame(Sample=mis[,1], Group=mis[,2])
mydf <- mis%>%
       group_by(Sample, Group)%>%
       summarize(Freq=n())
       
mydf <- mydf%>%
        group_by(Sample)%>%
        mutate(Perc=Freq/sum(Freq)*100, 
               Sample2=xSample2[as.character(Sample)])
               
fig0 <- ggplot(mydf,aes(x=as.character(Sample2), y=Group, fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("Individual")+ylab("Batch")+
            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
            
figfn <- "./1_demuxlet_output/figure2.mis.newIndi.pdf"
pdf(figfn, width=15, height=5)
fig0
dev.off()

###
### match data
match0 <- sample_batch[(sample_batch%in%correctBatch)]       
match0 <- str_split(match0, "_", simplify=TRUE)
match0 <- data.frame(Sample=match0[,1], Group=match0[,2])
mydf <- match0%>%
       group_by(Sample, Group)%>%
       summarize(Freq=n())
       
mydf <- mydf%>%
        group_by(Sample)%>%
        mutate(Perc=Freq/sum(Freq)*100, 
               Sample2=xSample2[as.character(Sample)])
               
fig0 <- ggplot(mydf,aes(x=as.character(Sample2), y=Group, fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("Individual")+ylab("Batch")+
            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
            
figfn <- "./1_demuxlet_output/figure3.match.indi.pdf"
pdf(figfn, width=15, height=5)
fig0
dev.off()


###
### all data
allx <- sample_batch      
allx <- str_split(allx, "_", simplify=TRUE)
allx <- data.frame(Sample=allx[,1], Group=allx[,2])
mydf <- allx%>%
       group_by(Sample, Group)%>%
       summarize(Freq=n())
       
mydf <- mydf%>%
        group_by(Sample)%>%
        mutate(Perc=Freq/sum(Freq)*100, 
               Sample2=xSample2[as.character(Sample)])
               
fig0 <- ggplot(mydf,aes(x=as.character(Sample2), y=Group, fill=Perc))+
            geom_tile()+
            scale_fill_gradient(low="white",high="#132B43",na.value=NA)+
            xlab("Individual")+ylab("Batch")+
            theme_bw()+theme(axis.text.x=element_text(angle=60,hjust=1))
            
figfn <- "./1_demuxlet_output/figure4.all.newIndi.pdf"
pdf(figfn, width=15, height=5)
fig0
dev.off()


               
     



