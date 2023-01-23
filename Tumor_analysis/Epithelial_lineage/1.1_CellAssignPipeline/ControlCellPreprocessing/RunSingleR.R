#SingleR Run


#Step 1: Convert to Seurat
library(Seurat)
library(SeuratData)
library(SeuratDisk)
Sys.setenv(RETICULATE_PYTHON = "~/Projects/Summer22/summer_venv/bin/python")
library(reticulate)
library(celldex)
library(SingleR)
use_virtualenv("~/Projects/Summer22/summer_venv/")
setwd("~/Projects/Summer22/Preprocessing/data/SampleQCCheckpoints")

ad <- import("anndata", convert = FALSE)
sc <- import("scanpy")


atlas.data <- sc$read_h5ad("ctrl_c51.h5ad")
counts <- t(atlas.data$X)
colnames(counts) <-  atlas.data$obs_names$to_list()
rownames(counts) <-  atlas.data$var_names$to_list()
counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
seurat <- CreateSeuratObject(counts)
seurat <- AddMetaData(seurat,  atlas.data$obs)
hpca.se <- BlueprintEncodeData()
hpca.se
pred.hesc <- SingleR(test = seurat, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)
c51_results <- SingleR(GetAssayData(seurat, assay = "RNA", slot = "data"), clusters = seurat@meta.data$leiden,ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.main)
write.csv(c51_results, "SingleRResults/c51_clusters.csv")


hpca.se <- BlueprintEncodeData()
hpca.se
pred.hesc <- SingleR(test = seurat, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)

samples = list(2,3,4,5,7)
for (x in samples) {
  atlas.data <- sc$read_h5ad(paste("ctrl_c5",toString(x),".h5ad", sep=""))
  counts <- t(atlas.data$X)
  colnames(counts) <-  atlas.data$obs_names$to_list()
  rownames(counts) <-  atlas.data$var_names$to_list()
  counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
  seurat <- CreateSeuratObject(counts)
  seurat <- AddMetaData(seurat,  atlas.data$obs)
  temp_results <- SingleR(GetAssayData(seurat, assay = "RNA", slot = "data"), clusters = seurat@meta.data$leiden,ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.main)
  write.csv(temp_results, paste("SingleRResults/c5", toString(x), "_clusters.csv",sep = ""))
  
  print(x)
}

