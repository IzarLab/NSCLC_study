# Author: Somnath Tagore, Ph.D. Title: Sample integration
# Script Name: integration_samples.R
# Last Updated: 01/01/2021

#!/usr/bin/env Rscript


library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(purrr)
library(cowplot)
#BiocManager::install("DropletUtils")
library(DropletUtils)
#BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("SeuratObject")
#library(SeuratObject)

library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(Matrix)
#BiocManager::install("celldex")
library(celldex)

pat <- 'Integrated'

PA001 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA001_seurat.rds")
DefaultAssay(PA001) <- 'RNA'
PA001$orig.ident<-rep('PA001',length(PA001$orig.ident))

PA004 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA004_seurat.rds")
DefaultAssay(PA004) <- 'RNA'
PA004$orig.ident<-rep('PA004',length(PA004$orig.ident))

PA005 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA005_seurat.rds")
DefaultAssay(PA005) <- 'RNA'
PA005$orig.ident<-rep('PA005',length(PA005$orig.ident))

PA019 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA019_seurat.rds")
DefaultAssay(PA019) <- 'RNA'
PA019$orig.ident<-rep('PA019',length(PA019$orig.ident))

PA025 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA025_seurat.rds")
DefaultAssay(PA025) <- 'RNA'
PA025$orig.ident<-rep('PA025',length(PA025$orig.ident))

PA034 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA034_seurat.rds")
DefaultAssay(PA034) <- 'RNA'
PA034$orig.ident<-rep('PA034',length(PA034$orig.ident))

PA042 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA042_seurat.rds")
DefaultAssay(PA042) <- 'RNA'
PA042$orig.ident<-rep('PA042',length(PA042$orig.ident))

PA043 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA043_seurat.rds")
DefaultAssay(PA043) <- 'RNA'
PA043$orig.ident<-rep('PA043',length(PA043$orig.ident))

PA048 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA048_seurat.rds")
DefaultAssay(PA048) <- 'RNA'
PA048$orig.ident<-rep('PA048',length(PA048$orig.ident))

PA054 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA054_seurat.rds")
DefaultAssay(PA054) <- 'RNA'
PA054$orig.ident<-rep('PA054',length(PA054$orig.ident))

PA056 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA056_seurat.rds")
DefaultAssay(PA056) <- 'RNA'
PA056$orig.ident<-rep('PA056',length(PA056$orig.ident))

PA060 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA060_seurat.rds")
DefaultAssay(PA060) <- 'RNA'
PA060$orig.ident<-rep('PA060',length(PA060$orig.ident))

PA067 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA067_seurat.rds")
DefaultAssay(PA067) <- 'RNA'
PA067$orig.ident<-rep('PA067',length(PA067$orig.ident))

PA068 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA068_seurat.rds")
DefaultAssay(PA068) <- 'RNA'
PA068$orig.ident<-rep('PA068',length(PA068$orig.ident))

PA070 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA070_seurat.rds")
DefaultAssay(PA070) <- 'RNA'
PA070$orig.ident<-rep('PA070',length(PA070$orig.ident))

PA072 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA072_seurat.rds")
DefaultAssay(PA072) <- 'RNA'
PA072$orig.ident<-rep('PA072',length(PA072$orig.ident))

PA076 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA076_seurat.rds")
DefaultAssay(PA076) <- 'RNA'
PA076$orig.ident<-rep('PA076',length(PA076$orig.ident))

PA080 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA080_seurat.rds")
DefaultAssay(PA080) <- 'RNA'
PA080$orig.ident<-rep('PA080',length(PA080$orig.ident))

PA104 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA104_seurat.rds")
DefaultAssay(PA104) <- 'RNA'
PA104$orig.ident<-rep('PA104',length(PA104$orig.ident))

PA125 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA125_seurat.rds")
DefaultAssay(PA125) <- 'RNA'
PA125$orig.ident<-rep('PA125',length(PA125$orig.ident))

PA141 <- readRDS(file="/home/ubuntu/NSCLC/seurat_objs/PA141_seurat.rds")
DefaultAssay(PA141) <- 'RNA'
PA141$orig.ident<-rep('PA141',length(PA141$orig.ident))

nsclc.combined <- merge(PA001, y =c(PA004,PA005,PA019,PA025,PA034,PA042,PA043,PA048,PA054,PA056,PA060,PA067,PA068,PA070,PA072,PA076,PA080,PA104,PA125,PA141), add.cell.ids = c('PA001','PA004','PA005','PA019','PA025','PA034','PA042','PA043','PA048','PA054','PA056','PA060','PA067','PA068','PA070','PA072','PA076','PA080','PA104','PA125','PA141'), project = 'NSCLC')

#nsclc.combined <- merge(PA001, y =c(PA004,PA005), add.cell.ids = c('PA001','PA004','PA005'), project = 'NSCLC')

nsclc.combined

table(nsclc.combined$orig.ident)

nsclc.combined.list <- SplitObject(nsclc.combined, split.by = "orig.ident")
for (i in 1:length(nsclc.combined.list)) {
  print(i)
  #nsclc.combined.list[[i]] <- NormalizeData(nsclc.combined.list[[i]], verbose = FALSE , scale.factor = 1000000)
  nsclc.combined.list[[i]] <- NormalizeData(nsclc.combined.list[[i]],verbose=FALSE)
  #nsclc.combined.list[[i]] <- FindVariableFeatures(nsclc.combined.list[[i]], selection.method = "vst", nfeatures = 2000)
  nsclc.combined.list[[i]] <- FindVariableFeatures(nsclc.combined.list[[i]], verbose=FALSE)
}
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = nsclc.combined.list)
nsclc.combined.list <- lapply(X = nsclc.combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find integration anchors
#anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30, reduction = "rpca")

# Integrate data sets
#seu <- IntegrateData(anchorset = anchors, dims = 1:30)
#saveRDS(seu, file = paste0("/home/ubuntu/SCLC/",seu,'_integrated_samples_rpca.rds'))

#features <- SelectIntegrationFeatures(object.list = sclc.combined.list,)
#patients.anchors <- FindIntegrationAnchors(object.list = nsclc.combined.list, anchor.features = features)
patients.anchors <- FindIntegrationAnchors(object.list = nsclc.combined.list,  reduction = "rpca",dims = 1:30)
#patients.anchors <- FindIntegrationAnchors(object.list = sclc.combined.list, dims = 1:30, anchor.features = features)

#patients.anchors <- FindIntegrationAnchors(object.list = nsclc.combined.list, dims = 1:30, anchor.features = 2000)
patients.integrated <- IntegrateData(anchorset = patients.anchors,dims=1:30)
#patients.integrated <- IntegrateData(anchorset = patients.anchors)
#saveRDS(patients.integrated, file = paste0("/home/ubuntu/SCLC/",pat,'_integrated_after_CB_QC_seurat.rds'))
#saveRDS(patients.integrated, file = paste0("/home/ubuntu/NSCLC/seurat_objs/",pat,'_NSCLC_samples_after_QC_cb_seurat.rds'))
saveRDS(patients.integrated, file = paste0("/home/ubuntu/NSCLC/seurat_objs/",pat,'_NSCLC_samples_after_QC_cb_seurat_refined.rds'))

     
