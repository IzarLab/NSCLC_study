# Author: Somnath Tagore, Ph.D. Title: T-cells subclustering
# Script Name: tcells_clustering.R 
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

pat <- 'Integrated'
group <- 'tcells'
res <- 0.2

# Read the data
patients.integrated<-readRDS(file= paste0("Integrated_nontumor_44_samples_",group,"_after_manual_annotation_update_1.rds"))

pdf(file = paste0(,pat,"Integrated_nontumor_44_samples_",group,"_after_manual_annotation_update_1_v3.pdf"),height = 10,width=10)
#
print(dim(patients.integrated@assays$RNA@counts))
print('Number of genes')
print(dim(patients.integrated@assays$RNA@counts)[1])
print('Number of cells')
print(dim(patients.integrated@assays$RNA@counts)[2])
print('Median number of genes')
print(round(median(patients.integrated@meta.data$nFeature_RNA)))

immune.combined<-patients.integrated

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = res)

cluster_0 <- subset(immune.combined,idents=0)
cluster_1 <- subset(immune.combined,idents=1)
cluster_2 <- subset(immune.combined,idents=2)
cluster_3 <- subset(immune.combined,idents=3)
cluster_4 <- subset(immune.combined,idents=4)
cluster_5 <- subset(immune.combined,idents=5)
cluster_6 <- subset(immune.combined,idents=6)
cluster_7 <- subset(immune.combined,idents=7)
cluster_8 <- subset(immune.combined,idents=8)
saveRDS(cluster_0,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_0_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_1,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_1_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_2,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_2_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_3,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_3_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_4,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_4_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_5,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_5_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_6,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_6_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_7,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_7_after_manual_annotation_update_1_v3.rds"))
saveRDS(cluster_8,paste0("/",group,"/",pat,"Integrated_nontumor_44_samples_",group,"_cluster_8_after_manual_annotation_update_1_v3.rds"))

# Visualization
DimPlot(object = immune.combined, reduction = "pca",label=T, label.size=7.2)+labs(title = paste0("PCA All 44 Samples Non-tumor cells: NSCLC_KRAS_STK-",pat,"\n",group," Resolution parameter =", res))
DimPlot(object = immune.combined, reduction = "umap",label=T, label.size=7.2)+labs(title = paste0("UMAP All 44 Samples Non-tumor cells: NSCLC_KRAS_STK-",pat,"\n",group," Resolution parameter =", res))

dev.off()
print(paste('End:',Sys.time()))

