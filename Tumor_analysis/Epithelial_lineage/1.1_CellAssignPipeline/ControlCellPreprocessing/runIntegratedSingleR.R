#Step 1: Convert to Seurat
library(Seurat)
library(SeuratData)
library(SeuratDisk)
Sys.setenv(RETICULATE_PYTHON = "~/Projects/Summer22/summer_venv/bin/python")
library(reticulate)
library(celldex)
library(SingleR)
library(sctransform)

use_virtualenv("~/Projects/Summer22/summer_venv/")
setwd("~/Projects/Summer22/Preprocessing/")



ad <- import("anndata", convert = FALSE)
sc <- import("scanpy")


atlas.data <- sc$read_h5ad("data/IntegratedDataQCCheckpoints/integrated_normalized_control.h5ad")
counts <- t(atlas.data$X)
colnames(counts) <-  atlas.data$obs_names$to_list()
rownames(counts) <-  atlas.data$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
seurat <- CreateSeuratObject(counts, meta.data = meta_dat)
seurat <- AddMetaData(seurat,  atlas.data$obs)
meta_dat <- atlas.data$obs




hpca.se <- BlueprintEncodeData()
hpca.se
temp_results <- SingleR(GetAssayData(seurat, assay = "RNA", slot = "data"), clusters = seurat@meta.data$leiden,ref = hpca.se, assay.type.test=1,
                        labels = hpca.se$label.main)

write.csv(temp_results, paste("data/IntegratedDataQCCheckpoints/integratedSingleR", "Clusters.csv",sep = ""))

