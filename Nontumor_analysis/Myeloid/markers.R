# Author: Somnath Tagore, Ph.D. Title: Markers for Myeloid cells
# Script Name: markers.R 
# Last Updated: 05/20/2021

# Packages required for this analysis
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(purrr)
library(cowplot)
library(DropletUtils)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(Matrix)
library(celldex)
library(matrixStats)
library(copynumber)
library(infercnv)
library(stringr)


library(tidyverse)
library(hrbrthemes)
library(viridis)


library("ggpubr")
library("ggExtra")

library(cowplot)


pat <- 'Group3_Myeloid_cells'

# Color palette v2
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')

seu_all <- readRDS('./myeloid/NSCLC_44_samples_Myeloid_v7.rds')
seu <- seu_all


DefaultAssay(seu)<-'RNA'

genes_myeloid_v10 <- c('CLEC9A',
                       'CLNK',
                       'ENPP1',
                       'CPVL',
                       'PAX7',
                       'MUC1',
                       'SFTPB',
                       'CDCA8',
                       'CENPF',
                       'ITGAM',
                       'CD274',
                       'LAMP3',
                       'SAMSN1',
                       'CD200',
                       'ITGAX',
                       'CAMK1D',
                       'ABCC5',
                       'CSF2RA',
                       'NAMPT',
                       'CD63',
                       'CD2',
                       'ITK',
                       'IKZF3',
                       'CD1A',
                       'CD207',
                       'CLDN1',
                       'HIP1',
                       'PARM1',
                       'STAT1',
                       'GBP1',
                       'GBP5',
                       'MX1',
                       'MX2',
                       'IFI44',
                       'SIGLEC1',
                       'RNASE1',
                       'CD163',
                       'MERTK',
                       'MARCO',
                       'CYP27A1',
                       'PPARG',
                       'ANO5',
                       'LTA4H',
                       'FBP1',
                       'DENND4C',
                       'ABCG1',
                       'PBX1',
                       'GATA2',
                       'KIT',
                       'IL18R1',
                       'MS4A2',
                       'P2RY12',
                       'NAV3',
                       'SLC2A5',
                       'NCK2',
                       'SPP1',
                       'VCAN',
                       'FCN1',
                       'VIM',
                       'LYZ',
                       'CD14')


pdf('./Myeloid_markers_dot.pdf', height = 15, width = 15)
# replace <- Vectorize(function(x) {
#   gsub("_"," ",x)
# })
replace <- Vectorize(function(x) {
  gsub("_"," ",x)
})
DotPlot(seu, features = unique(genes_myeloid_v10), group.by=  "manual.annot.fine", dot.scale = 5) +
#  DotPlot(seu, features = unique(genes_t_cells_v9), group.by=  "manual.annot.fine", dot.scale = 5) +
  RotatedAxis() + labs(y="Class", x = "Genes") + theme(axis.text=element_text(size=13)) + 
  scale_y_discrete(labels = replace) +
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90),
        panel.grid.major = element_line(colour = "gray",linetype = "dashed", size=0.35), 
        panel.border = element_rect(colour = "black", fill=NA, size=2))+coord_flip()+RotatedAxis()+
  facet_grid( ~ unique(seu@meta.data$PRIMARY_vs_BRAIN_METS))


dev.off()

pdf('./Myeloid_markers_violin.pdf', height = 10, width = 6)

VlnPlot(seu, features = genes_myeloid_v10, group.by = 'manual.annot.fine', pt.size = 0,
                assay = 'RNA', stack = T, flip = T,split.by = 'Primary_vs_BM',cols = colPrimvsBMmain) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom',
        axis.title.x = element_blank())


dev.off()

