# Author: Somnath Tagore, Ph.D. Title: Slideseq Reference creation
# Script Name: reference_creation.R 
# Last Updated: 10/01/2022

library(spacexr)
library(Matrix)
library(stringr)
library(Seurat)

datadir = "/home/ubuntu/NSCLC/slide-seq/"

#load in single-nuclei reference, filter for only single nuclei data, and filter out low-quality cells
data_nsclc = readRDS(paste0(datadir,"/sndata_all_samples.rds"))
DefaultAssay(data_nsclc) = "RNA"
# data_MBPM$placeholder = !(data_MBPM$cell_type_int %in% c("Low-quality cells","Doublets","Contamination","Undetermined"))
# data_MBPM = subset(data_MBPM, placeholder)
# data_MBPM = subset(data_MBPM, sequencing=="Single nuclei")

counts = data_nsclc@assays$RNA@counts
genenames = rownames(data_nsclc)
barcodes = colnames(data_nsclc)
rownames(counts) = genenames
colnames(counts) = barcodes

#extract cell_type_main information from single nuclei data object, create Reference object for rctd
cell_types = data_nsclc$tumor_nontumor_categories
#cell_types = data_nsclc$cell_type_int
data_nsclc$cell_type_broad <- data_nsclc$tumor_nontumor_categories
#data_MBPM$cell_type_broad = "Non-immune"
#data_MBPM$cell_type_broad[data_MBPM$cell_type_int %in% c("B cells","CD4+ T cells","CD8+ T cells","Dendritic cells","Mast cells","MDM","Monocytes","NK cells","Plasma cells","Tregs")] = "Immune"
#data_MBPM$cell_type_broad[data_MBPM$cell_type_int %in% c("Tumor cells")] = "Tumor"
Idents(data_nsclc) = data_nsclc$cell_type_broad
#cell_types = data_MBPM$cell_type_broad
#cell_types = str_replace_all(cell_types,"/","_")
names(cell_types) = barcodes
cell_types = factor(cell_types)

reference <- Reference(counts, cell_types)
saveRDS(file="reference.slideseq.rds")
