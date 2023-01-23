# Author: Somnath Tagore, Ph.D. Title: Downstream analysis (Sample-wise) 
# Script Name: sample_wise.R
# Last Updated: 01/01/2021

#!/usr/bin/env Rscript

### Seurat analysis for cellbender output provided as an argument
## needs to run in /home/ubuntu

print(paste('Start:',Sys.time()))

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

#Input argument
pat <- 'STK_15'

doublet_rate <-0.0483
print(pat)
#print(doublet_rate)

##### Loading, merging, QC, dimension reduction #####
### Load dataset
seu.data <- Read10X_h5(paste0(pat,'_filtered.h5'), use.names = TRUE, unique.features = TRUE)

### Initialize the Seurat object with the raw (non-normalized data)
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 3, min.features = 200)

# Annotate MT genes
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^MT-", col.name = "percent.mt")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPS", col.name = "percent.rps")
seu_raw <- PercentageFeatureSet(seu_raw, pattern = "^RPL", col.name = "percent.rpl")
seu_raw$percent.rp=seu_raw$percent.rps + seu_raw$percent.rpl

# Identify doublets using scrublet
writeMM(seu_raw@assays$RNA@counts, paste0(pat,'_raw.mtx'))
#system(paste('/home/ubuntu/anaconda3/bin/python /home/ubuntu/SCLC/scrublet.run.py', pat))#, doublet_rate))
doublets <- read.table(paste0(pat,'_raw.txt'),header = T)
seu_raw[['predicted_doublets']]<-doublets$predicted_doublets
seu_raw[['doublet_scores']]<-doublets$doublet_scores
#system(paste0('rm data/',pat,'/matrix_',pat,'_raw.mtx'))
#system(paste0('rm data/',pat,'/doublets_',pat,'_raw.txt'))

### subset
minFeature<-200
maxFeature<- 9000
minCount<- 800
maxCount<- 40000
maxMT<-15
seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
                nCount_RNA < maxCount & percent.mt < maxMT & predicted_doublets ==F)
#seu <- subset(seu_raw, subset = percent.mt < maxMT & predicted_doublets ==F)

seu

### Workflow RNA
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

### Workflow SCT
seu <- SCTransform(seu,assay='RNA',vars.to.regress = c("nCount_RNA","percent.mt","percent.rp"))
seu <- RunPCA(seu,assay = 'SCT')
seu <- RunUMAP(seu, dims = 1:30, assay = 'SCT')
seu <- FindNeighbors(seu, dims = 1:30,assay='SCT')
seu<- FindClusters(seu)

### Finding differentially expressed features (cluster biomarkers)
#seu.markers <- FindAllMarkers(seu, only.pos = F, min.pct = 0.2, logfc.threshold = 0.25)
#marker_table <- seu.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#write.csv(marker_table,paste0("data/",pat,'/markers_',pat,'_cellbender.csv'),row.names = F)

### cell type identification
#seu_sce <- as.SingleCellExperiment(seu)
seu_sce<-seu
bped<-BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce@assays$RNA@counts, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']]<-pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce@assays$RNA@counts, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']]<-pred_bped_fine$pruned.labels

# iced<-DatabaseImmuneCellExpressionData()
# pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
# pruneScores(pred_iced_main)
# seu[['celltype_iced_main']]<-pred_iced_main$pruned.labels
# pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
# pruneScores(pred_iced_fine)
# seu[['celltype_iced_fine']]<-pred_iced_fine$pruned.labels

hpca<-HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce@assays$RNA@counts, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']]<-pred_hpca_main$pruned.labels
pred_hpca_fine <- SingleR(test = seu_sce@assays$RNA@counts, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']]<-pred_hpca_fine$pruned.labels


### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_raw_features','n_raw_cells','n_predicted_doublets','n_features','n_cells','median_features','median_counts','cutoff_features','cutoff_counts','cutoff_mt')
rownames(stats)<-pat
stats$sample<-pat
stats$n_raw_features<-dim(seu_raw@assays$RNA@counts)[1]
stats$n_raw_cells<-dim(seu_raw@assays$RNA@counts)[2]
stats$n_predicted_doublets <-length(which(seu_raw@meta.data$predicted_doublets ==T))
stats$n_features<-dim(seu@assays$RNA@counts)[1]
stats$n_cells<-dim(seu@assays$RNA@counts)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))
stats$cutoff_features<-paste(minFeature,maxFeature)
stats$cutoff_counts<-paste(minCount,maxCount)
stats$cutoff_mt<-paste(maxMT)



print(dim(seu@assays$RNA@counts)[1])
print(dim(seu@assays$RNA@counts)[2])
print(round(median(seu@meta.data$nFeature_RNA)))
saveRDS(seu, file = paste0("~/Documents/Ben_Izar_project/NSCLC/Seurat/v4.0_new_version/",pat,"/NSCLC_",pat,'_cb_after_QC_seurat.rds'))


### write pdf reports
pdf(file = paste0(pat,"_celltypes.pdf"),height=10,width=10)

# stats
textplot(t(stats),cex=1.2,halign='left')

# plots raw data
print(qplot(x=seu_raw@meta.data$nCount_RNA,y = seu_raw@meta.data$nFeature_RNA, col=seu_raw@meta.data$percent.mt, xlab = "nCount_RNA",
            ylab = "nFeature_RNA", main =paste0(seu_raw@meta.data$patient[1], " raw data: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = "Percent MT") + theme_classic())
print(VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,split.by=NULL))

# QC plot filtered
print(qplot(x=seu@meta.data$nCount_RNA,y = seu@meta.data$nFeature_RNA, col=seu@meta.data$percent.mt, xlab = "nCount_RNA",
            ylab = "nFeature_RNA", main =paste0(seu@meta.data$patient[1], " filtered: nCount_RNA vs. nFeature_RNA")) + scale_colour_gradient(low="blue", high="green") + labs(color = 'Percent MT') + theme_classic())
print(VlnPlot(seu, features = 'nFeature_RNA'))
print(VlnPlot(seu, features = "nCount_RNA"))
print(VlnPlot(seu, features = "percent.mt"))

# Determine metrics to plot present in seurat_control@meta.data
metrics <-  c("nCount_SCT", "nFeature_SCT", "percent.mt",'doublet_scores')
# Extract the UMAP coordinates for each cell and include information about the metrics to plot
qc_data <- FetchData(seu, vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seu, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  as.data.frame() %>% group_by(ident) %>% summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plot a UMAP plot for each metric
map(metrics, function(qc){
  ggplot(qc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=qc), size=0.5) +
    scale_color_gradientn(colours = rev(rainbow(4))) +theme_classic() +
    geom_text(data=umap_label, aes(label=ident, x, y)) +
    ggtitle(qc)+
    theme(legend.text=element_text(size=10),legend.title=element_text(size=14),legend.key.width=unit(0.2,"cm"))
}) %>% plot_grid(plotlist = .)


# PCA
print(DimPlot(seu, reduction = "pca",group.by = 'ident',label.size=5.2))
print(DimPlot(seu, reduction = "pca",group.by = 'Cell_types_manual',label = T,label.size=5.2, cols = c("Tumor"="red","B_cells"="forestgreen", "Fibroblasts"="coral4","Myeloid"="cyan","Endothelial"="chartreuse4","T_cells"="burlywood4","CNS"="goldenrod1","De-differentiated_cells"="gray"))+ggtitle('STK_15: Celltypes'))

# UMAP
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',label.size=5.2))
print(DimPlot(seu, reduction = "umap",group.by = 'Cell_types_manual',label = T,label.size=5.2,cols = c("Tumor"="red","B_cells"="forestgreen", "Fibroblasts"="coral4","Myeloid"="cyan","Endothelial"="chartreuse4","T_cells"="burlywood4","CNS"="goldenrod1","De-differentiated_cells"="gray"))+ggtitle('STK_15: Celltypes'))


## bped
plotScoreHeatmap(pred_bped_fine, clusters=seu_sce@colData$ident,fontsize = 10,main='pred_bped_fine')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,label.size = 5.2) +
  ggtitle('celltype_bped_fine') +
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=10))

## hpca
plotScoreHeatmap(pred_hpca_main, clusters=seu_sce@colData$ident,fontsize = 10,main='pred_hpca_main')
DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,label.size = 5.2) +
  ggtitle('celltype_hpca_main') +
  guides(col = guide_legend(nrow = 30,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=10))

dev.off()


print(paste('End:',Sys.time()))
