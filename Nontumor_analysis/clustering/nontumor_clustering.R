# Author: Somnath Tagore, Ph.D. Title: Non-tumor subclustering for each major groups 
# Script Name: nontumor_clustering.R 
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


pat <- 'Group_1'

# Read the integerated data
patients.integrated<-readRDS(file= "Integrated_NSCLC_KRAS_STK_samples_nontumor_44_samples.rds")

DefaultAssay(patients.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
patients.integrated <- ScaleData(patients.integrated, verbose = FALSE)
patients.integrated <- RunPCA(patients.integrated, npcs = 30, verbose = FALSE)
patients.integrated <- RunUMAP(patients.integrated, reduction = "pca", dims = 1:30)
patients.integrated <- FindNeighbors(patients.integrated, reduction = "pca", dims = 1:30)
patients.integrated <- FindClusters(patients.integrated, resolution = 0.07)

pdf(file = paste0(pat,"_clustering.pdf"),height = 10,width=10)

# Group_1: Create PCA/UMAP
DimPlot(object = patients.integrated, reduction = "pca",label=T, label.size=7.2)+labs(title = paste0("PCA Non-tumor cells,",pat))
DimPlot(object = patients.integrated, reduction = "umap",label=T, label.size=7.2)+labs(title = paste0("UMAP Non-tumor cells,",pat))

textplot('Group_1',cex=1.2,halign = "center", valign = "center")
seu <- patients.integrated

# Plot relevant markers
if("CDH5" %in% rownames(seu@assays$RNA@data) |
   "CLDN5" %in% rownames(seu@assays$RNA@data) |
   "CLEC14A" %in% rownames(seu@assays$RNA@data)|
   "CXorf36" %in% rownames(seu@assays$RNA@data)|
   "ECSCR" %in% rownames(seu@assays$RNA@data)|
   "F2RL3" %in% rownames(seu@assays$RNA@data)|
   "FLT1" %in% rownames(seu@assays$RNA@data)|
   "FLT4" %in% rownames(seu@assays$RNA@data)|
   "GPR4" %in% rownames(seu@assays$RNA@data)|
   "GPR182" %in% rownames(seu@assays$RNA@data)|
   "KDR" %in% rownames(seu@assays$RNA@data)|
   "MMRN1" %in% rownames(seu@assays$RNA@data)|
   "MMRN2" %in% rownames(seu@assays$RNA@data)|
   "MYCT1" %in% rownames(seu@assays$RNA@data)|
   "PTPRB" %in% rownames(seu@assays$RNA@data)|
   "RHOJ" %in% rownames(seu@assays$RNA@data)|
   "SLCO2A1" %in% rownames(seu@assays$RNA@data)|
   "SOX18" %in% rownames(seu@assays$RNA@data)|
   "STAB2" %in% rownames(seu@assays$RNA@data)|
   "VWF" %in% rownames(seu@assays$RNA@data)   ) {
  print(FeaturePlot(seu, features = c("CDH5","CLDN5","CLEC14A","CXorf36", "ECSCR","F2RL3","FLT1", "FLT4","GPR4","GPR182","KDR","MMRN1" ,"MMRN2","MYCT1","PTPRB","RHOJ","SLCO2A1","SOX18","STAB2","VWF"), min.cutoff = "q05",max.cutoff = 'q95',order=T))}


dev.off()

