###
library("rhdf5")
library("corpcor")
library(Matrix)
library(MASS)
library(scales)
library(tidyverse) 
library(parallel)
library(data.table)
library(future)
library(purrr) 
library(furrr) 
library("BiocParallel")  
library(Rcpp)  
library(reshape)
library(qqman)   
library(qvalue)
##
library(Seurat)
library(SeuratDisk)
library(harmony)
library(cellranger)
library(diem)
library(DoubletFinder) 
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
library(ggsci)
library(gtable)
library(gg3D)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(UpSetR)
library(ComplexHeatmap)
library(gtable)
library(RColorBrewer) 
library(viridis)       
theme_set(theme_grey())

##constant
#MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
#Contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")

#col1 <- c("CTRL"="#828282", 
#           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
#           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")
#col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", 
#          "NKcell"="#aa4b56", "Tcell"="#ffaa00")
          
#lab2treat <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX", 
#               "PHA"="PHA", "PHA-DEX"="PHA+DEX")

          
                    
