# Author: Somnath Tagore, Ph.D. Title: Integrated Myeloid cells diffusion
# Script Name: myeloid_clustering.R 
# Last Updated: 05/20/2021

# Packages required for this analysis
library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)
library(scales)
#library(hypeR)


# Color palette
colors <- c('#006E82', '#AA0A3C', '#F0E442', '#00A0FA', '#FA5078', '#005AC8', '#CC79A7', '#FAE6BE', '#0072B2', '#A0FA82', '#F0F032', '#0AB45A',
            '#FA7850', '#14D2DC', '#FA78FA')

colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','Epithelial cells'='#FF00FF','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Epithelial cells'='#FF00FF', 'Non-STK11-mut'='#FFFF00', 'STK11-mut'='#000000')





# Read in file
patient <- readRDS('./Diffusion_Component_Analysis/Myeloid/scrnaseq_all_upd.rds')

colnames(patient@meta.data)

# Subset to .... cells

table(patient@meta.data$manual.annot)

####Cell cycle scoring

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps

patient <- FindVariableFeatures(patient, selection.method = "vst")
patient <- ScaleData(patient, features = rownames(patient))

patient <- RunPCA(patient)

patient <- CellCycleScoring(patient, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

head(patient[[]])

pdf('Myeloid_cells_res_0.1_refined_upd.pdf')
RidgePlot(patient, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

patient <- RunPCA(patient, features = c(s.genes, g2m.genes))
DimPlot(patient)

patient <- ScaleData(patient, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(patient))

patient <- RunPCA(patient, features = VariableFeatures(patient), nfeatures.print = 10)

patient <- RunPCA(patient, features = c(s.genes, g2m.genes))
DimPlot(patient)

dev.off()


#####

dim(patient@meta.data)

colnames(patient@meta.data)

# Subset to .... cells

table(patient@meta.data$cloneType)


DefaultAssay(patient) <- 'RNA'

# Rerun Seurat workflow
#patient <- ScaleData(object = patient)
#patient <- RunPCA(object = patient)
patient <- NormalizeData(patient)
patient <- FindVariableFeatures(patient, selection.method = "vst")
patient <- ScaleData(patient, features = rownames(patient))

patient <- RunPCA(patient)

patient <- RunUMAP(object = patient, dims = 1:30)

# Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(patient@assays$RNA@data)))

es@phenoData@data <- patient@meta.data

# Make diffusion map
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)
saveRDS(dm,file="Myeloid_cells_res_0.1_refined_upd_dm.rds")

# Add DCs to Seurat object
patient <- AddMetaData(patient, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

pdf('./Myeloid_cells_res_0.1_refined_upd.pdf')

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
#palette(c('#0072B2', '#0AB45A'))
plot(dm, col = as.factor(es@phenoData@data$STK11mut_vs_NonSTK11mut), main = 'STK11mut_vs_NonSTK11mut', pch = 16)
#legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$celltype_bped_main)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$celltype_bped_main))),
legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$STK11mut_vs_NonSTK11mut)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$STK11mut_vs_NonSTK11mut))),
       bty = 'n')

plot(dm, col = as.factor(es@phenoData@data$Primary_vs_BM), main = 'Primary_vs_BM', pch = 16)
#legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$celltype_bped_main)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$celltype_bped_main))),
legend('bottomright', inset = c(0, 0), legend = levels(as.factor(es@phenoData@data$Primary_vs_BM)), pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$Primary_vs_BM))),
       bty = 'n')


dev.off()


dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:3]

pdf('Myeloid_cells_res_0.1_refined_upd_2Dversion.pdf')

# Diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)

#palette(colPrimvsBMmain)


ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'manual.annot']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('manual.annot')+
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'seurat_clusters']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('seurat_clusters')+
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'celltype_bped_main']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('celltype_bped_main')+
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'celltype_bped_fine']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('celltype_bped_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'nCount_RNA']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('nCount_RNA')+
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'nFeature_RNA']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('nFeature_RNA')+
  theme(legend.title = element_blank())


palette(colPrimvsBMmain)
ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'Primary_vs_BM']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('Primary_vs_BM')+
  theme(legend.title = element_blank())


palette(colSTKvsNonSTKmain)
ggplot(dm_df, aes(DC1, DC2, col = patient@meta.data[rownames(dm_df), 'STK11mut_vs_NonSTK11mut']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('STK11mut_vs_NonSTK11mut')+
  theme(legend.title = element_blank())


dev.off()
