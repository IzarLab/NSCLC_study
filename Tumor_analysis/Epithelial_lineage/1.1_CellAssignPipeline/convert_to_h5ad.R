library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Matrix)


rds_temp <- readRDS("Integrated_NSCLC_KRAS_STK_samples_tumor_44_samples_annotations.rds")
matrix = GetAssayData(object = rds_temp, slot = "counts", assay="RNA")

writeMM(GetAssayData(object = rds_temp, slot = "counts", assay="RNA"), 'IntegratedTumorsRNARaw.mtx')
