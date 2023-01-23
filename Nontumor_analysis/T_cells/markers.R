# Author: Somnath Tagore, Ph.D. Title: Markers for T cells
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


pat <- 'Group1_T_cells'

# Color palette v2
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')

seu_all <- readRDS('./myeloid/NSCLC_44_samples_Tcells_v7.rds')
seu <- seu_all


DefaultAssay(seu)<-'RNA'

genes_t_cells_v9 <- c('BCL11B',
                      'FYN',
                      'IKZF1',
                      'CCL5',
                      'LYST',
                      'CD247',
                      'CBLB',
                      'IQGAP2',
                      'ETS1',
                      'SKAP1',
                      'ELMO1',
                      'ITGA4',
                      'DTDH1',
                      'AOAH',
                      'PTPRC',
                      'CD44',
                      'PYHIN1',
                      'CD8A',
                      'CD96',
                      'IFI6',
                      'NFAT5',
                      'THEMIS',
                      'PRKCH',
                      'ITK',
                      'PDE3B',
                      'BACH2',
                      'NCK2',
                      'CELF2',
                      'ITGA1',
                      'RBPJ',
                      'TNIK',
                      'LEF1',
                      'RIPOR2',
                      'AKT3',
                      'RORA',
                      'FOXP1',
                      'NCAM1',
                      'KLRC1',
                      'KLRC2',
                      'KLRD1',
                      'LYN',
                      'NCR1',
                      'KLRK1',
                      'RPLP1',
                      'RPL41',
                      'RPL13',
                      'MT-CO2',
                      'ACTB',
                      'B2M',
                      'STAT4',
                      'RUNX2',
                     'ITGA4',
                     # 'IKZF1',
                      'ITGA1',
                      'CD69',
                      'ITGAE',
                      'CD44',
                      'CD58',
                      'IL2RA',
                      'TBC1D4',
                      'CTLA4',
                      'ICOS',
                      'CD4')

pdf('./Tcells_markers_dot.pdf', height = 15, width = 15)
# replace <- Vectorize(function(x) {
#   gsub("_"," ",x)
# })
replace <- Vectorize(function(x) {
  gsub("_"," ",x)
})
DotPlot(seu, features = unique(genes_t_cells_v9), group.by=  "manual.annot.fine", dot.scale = 5) +
  RotatedAxis() + labs(y="Class", x = "Genes") + theme(axis.text=element_text(size=13)) + 
  scale_y_discrete(labels = replace) +
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90),
        panel.grid.major = element_line(colour = "gray",linetype = "dashed", size=0.35), 
        panel.border = element_rect(colour = "black", fill=NA, size=2))+coord_flip()+RotatedAxis()+
  facet_grid( ~ unique(seu@meta.data$PRIMARY_vs_BRAIN_METS))


dev.off()

pdf('./Tcells_markers_violin.pdf', height = 10, width = 6)

VlnPlot(seu, features = genes_t_cells_v9, group.by = 'manual.annot.fine', pt.size = 0,
                assay = 'RNA', stack = T, flip = T,split.by = 'Primary_vs_BM',cols = colPrimvsBMmain) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom',
        axis.title.x = element_blank())


dev.off()
