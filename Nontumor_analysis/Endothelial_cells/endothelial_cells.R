#!/usr/bin/env Rscript

### title: Analysis of endothelial cells from NSCLC
### author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(destiny)
library(SingleCellExperiment)
library(ggpubr)
'%notin%' <- Negate('%in%')

path.ct<-'data/endothelial/'
filename<-'NSCLC_endothelial'

colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colPrimvsBMmain2 <- c('#FF0000','yellow','#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')


setwd('~/Documents/NSCLC/')

system(paste0("aws s3 sync s3://nscl-seq/Seurat/integrated_data/44_samples/nontumor_jana/endothelial/ ~/Documents/NSCLC/data/endothelial/ --exclude '*' --include '*after_manual_annotation_v5*' --quiet"))

seu<- readRDS(paste0(path.ct,'Nontumor_group7_after_manual_annotation_v5.rds'))
Idents(seu)<-seu$integrated_snn_res.0.1
seu$barcode<-rownames(seu@meta.data)
seu <- FindClusters(seu, resolution = 0.1)


# get module scores for sigs
sigs<-read.csv('~/brain_mets/signatures/endothelial.csv',na.strings = '')
for(c in 1:ncol(sigs)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs[,c])),
                      name = colnames(sigs)[c],assay = 'RNA',search=F)
}

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')


# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA', min.pct = 0.25, 
                          logfc.threshold = 0.25,test.use = 'MAST')
write.csv(markers, paste0(path.ct,'markers_',filename,'.csv'), row.names = F)
markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC) -> tops

DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu)<-'integrated'

pdf(file = paste0(path.ct,'markers_',filename,'_heatmap.pdf'),width = 15, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'ident',raster = T,assay = 'RNA')
dev.off()


# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- as.character(seu$manual.annot)
seu$cell_type_fine <- ifelse(seu$cell_type_fine %in% c("Mixed bag"),'Capillary ECs', 
                             seu$cell_type_fine)

seu$cell_type_main <- 'Endothelial cells'


### stats
stats<-as.data.frame(matrix(data=NA,nrow = 1, ncol = 11))
colnames(stats)<-c('sample','n_features','n_cells','median_features','median_counts')
rownames(stats)<-filename
stats$sample<-filename
stats$n_features<-dim(seu@assays$integrated@data)[1]
stats$n_cells<-dim(seu@assays$integrated@data)[2]
stats$median_features<-round(median(seu@meta.data$nFeature_RNA))
stats$median_counts<-round(median(seu@meta.data$nCount_RNA))

pdf(paste0(path.ct,'plots_',filename,'.pdf'))
textplot(t(stats),cex=1.2,halign='left')

print(DimPlot(seu, reduction = "umap",label = T,group.by = 'ident',shuffle = T,raster = T))
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'orig.ident',shuffle = T,raster = T))
DimPlot(seu,reduction = 'pca',group.by = 'cell_type_fine',shuffle = T,raster = T)
DimPlot(seu,group.by = 'cell_type_fine',shuffle = T,raster = T)
DimPlot(seu,group.by = 'cell_type_main',shuffle = T,raster = T)
FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle', pt.size = 0.01)
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'cell_cycle',shuffle = T,raster = T))

DimPlot(seu, reduction = "umap",label = T,group.by = 'manual.annot',repel = T,
        label.size = 3,shuffle = T,raster = T)
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T,cols = colSTKvsNonSTKmain)
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T,cols = colPrimvsBMmain)

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 3,shuffle = T,raster = T)+ 
  guides(col = guide_legend(nrow = 25,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,shuffle = T,raster = T)+ 
  guides(col = guide_legend(nrow = 25,override.aes = list(size=5))) +
  theme(legend.text=element_text(size=6))

DotPlot(seu,assay = 'RNA',features = paste0(names(sigs),1),
        group.by = 'integrated_snn_res.0.1', dot.scale = 3)+ 
  scale_color_viridis() + coord_flip()

DotPlot(seu,assay = 'RNA',features = paste0(names(sigs),1),
        group.by = 'cell_type_fine', dot.scale = 3)+ 
  scale_color_viridis() + coord_flip() +RotatedAxis()

print(FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))
print(FeaturePlot(seu, features = c("percent.rps", "percent.rpl", 'G2M.Score', "S.Score"), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))

for(i in seq(1, ncol(sigs), 4)){
  four<-paste0(colnames(sigs)[i:(i+3)],'1')
  print(FeaturePlot(seu, features = c(four),min.cutoff = "q05", max.cutoff = "q95",
                    order=T,raster=T))
}

VlnPlot(seu,c("nCount_RNA", "nFeature_RNA", 'percent.mt', "doublet_scores","percent.rps", "percent.rpl", 'G2M.Score', "S.Score"),
        pt.size = 0)+NoLegend()
dev.off()


# Save object and labels
saveRDS(seu,paste0(path.ct,'data_',filename,'.rds'))
ct <- seu@meta.data %>% dplyr::select('barcode', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct,'data_',filename,'_celltype.csv'))

seu<-readRDS(paste0(path.ct,'data_',filename,'.rds'))

# Differential gene expression (DGE) based on type
Idents(seu)<-seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA', min.pct = 0.25, 
                          logfc.threshold = 0.25,test.use = 'MAST')
write.csv(markers, paste0(path.ct,'markers_',filename,'_celltype.csv'), row.names = F)
markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(15, wt=avg_log2FC) -> tops
DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu)<-'integrated'
pdf(file = paste0(path.ct,'markers_',filename,'_heatmap_celltype.pdf'),width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine',raster = T,assay = 'RNA')
dev.off()


#### Diffusion v2 #####
# Read-in object 
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
seu <- subset(seu, cell_type_fine %notin% c('Lymphatic ECs') & seu$DC3> (-0.05))
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Diffusion pseudotime calculation
dpt <- DPT(dm)
es$pseudotime_dpt <- rank(dpt$dpt)
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

saveRDS(dm,paste0(path.ct,'data_',filename,'_dm.rds'))

# Plots
pdf(paste0(path.ct,'plots_',filename,'_diffusion.pdf'))

# UMAPs
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T,cols = colSTKvsNonSTKmain)
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T,cols = colSTKvsNonSTKmain) + NoLegend()
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T,cols = colPrimvsBMmain)
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T,cols = colPrimvsBMmain) + NoLegend()
DimPlot(seu, group.by = 'orig.ident', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'orig.ident', raster = T, shuffle = T) + NoLegend()

# 3D diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSTKvsNonSTKmain)
tmp_var<-seu$STK11mut_vs_NonSTK11mut
plot(dm, col = as.factor(tmp_var), main = 'STK11mut_vs_NonSTK11mut', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colPrimvsBMmain2)
tmp_var<-seu$PRIMARY_vs_BRAIN_METS
plot(dm, col = as.factor(tmp_var), main = 'PRIMARY_vs_BRAIN_METS', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette('ggplot2')
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette('ggplot2')
tmp_var<-seu$orig.ident
plot(dm, col = as.factor(tmp_var), main = 'orig.ident', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<- seu$pseudotime_dpt
plot(dm, col = sig_col, main = 'pseudotime_dpt', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)),
            posx = c(0.88, 0.9), posy = c(0, 0.65))

sig_col<-seu@assays$RNA@data['VEGFC',]
plot(dm, col = sig_col, main = 'VEGFC', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

sig_col<-seu@assays$RNA@data['INSR',]
plot(dm, col = sig_col, main = 'INSR', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

sig_col<-seu@assays$RNA@data['IL1R1',]
plot(dm, col = sig_col, main = 'IL1R1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

# 2D diffusion maps
set.seed(1)
dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))), 1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('cell_type_fine') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'STK11mut_vs_NonSTK11mut'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('STK11mut_vs_NonSTK11mut') + 
  theme(legend.title = element_blank())+
  scale_color_manual(values = colSTKvsNonSTKmain)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'PRIMARY_vs_BRAIN_METS'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('PRIMARY_vs_BRAIN_METS') + 
  theme(legend.title = element_blank())+
  scale_color_manual(values = colPrimvsBMmain2)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'orig.ident'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('orig.ident') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'pseudotime_dpt'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('pseudotime_dpt') + 
  theme(legend.title = element_blank()) + scale_color_viridis()

FeaturePlot(seu,features = paste0('DC',seq(1,9)),order = T,raster = F,pt.size = 0.01) & 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))

dev.off()

#### Heatmap ####
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
dm<-readRDS(paste0(path.ct,'data_',filename,'_dm.rds'))
seu <- subset(seu, barcode %in% rownames(dm@eigenvectors))
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))
DefaultAssay(seu)<-'RNA'
seu<-ScaleData(seu)

markers<-read.csv(paste0(path.ct,'markers_',filename,'.csv'))
markers %>% filter(p_val_adj < 0.05 & cluster %in% c(3,0,1) ) %>% 
  group_by(cluster) %>% top_n(50, wt=avg_log2FC)  -> tops
tops$cluster<-factor(tops$cluster,levels = c(3,0,1))
tops<-tops %>% arrange(cluster)
sel_genes<-tops$gene

seu_sel<-subset(seu,features = sel_genes)

#datamat<- rescale(seu_sel@assays$RNA@scale.data, to=c(-1,1))
datamat<- seu_sel@assays$RNA@scale.data
datamat[datamat>2]<-2
datamat = datamat[match(sel_genes, rownames(datamat)), ]
datamat = datamat[, order(seu_sel$pseudotime_dpt,decreasing = F)]

# marker anno
annotation_row = data.frame(DEG = c(rep('Arterial DEGs',50), rep('Capillary DEGs',50),
                                    rep('Venous DEGs',50)))
rownames(annotation_row) = rownames(datamat)


# cluster anno
cell_type_fine = seu_sel$cell_type_fine
cell_type_fine = cell_type_fine[order(seu_sel$pseudotime_dpt,decreasing = F)]
STK11vNON = seu_sel$STK11mut_vs_NonSTK11mut
STK11vNON = STK11vNON[order(seu_sel$pseudotime_dpt,decreasing = F)]
PRIvBM = seu_sel$PRIMARY_vs_BRAIN_METS
PRIvBM = PRIvBM[order(seu_sel$pseudotime_dpt,decreasing = F)]
pseudo = seu_sel$pseudotime_dpt
pseudo = pseudo[order(seu_sel$pseudotime_dpt,decreasing = F)]
pseudo<- rescale(pseudo, to=c(0,1))

annotation_col = data.frame(cell_type_fine = factor(cell_type_fine),
                            STK11vNON=factor(STK11vNON),
                            PRIvBM=factor(PRIvBM),
                            DPT=pseudo)
rownames(annotation_col) = colnames(datamat)

ann_colors = list(
  cell_type_fine = c(`Capillary ECs`="#F8766D",`Arterial ECs`="#00BFC4",
                     `Venous ECs`="#7CAE00",`Cycling ECs`="#C77CFF"),
  STK11vNON =c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000'),
  PRIvBM=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF','CHEST_WALL_MET'='yellow'),
  DPT=viridis(100))
 
hue_pal()(4)

pdf('data/endothelial/plots_heatmap_pseudotime_top50.pdf',height = 18,width = 10)
pheatmap(datamat,cluster_rows = F, cluster_cols = F, color = PurpleAndYellow(),
         show_colnames = F,annotation_col = annotation_col,annotation_row = annotation_row,
         annotation_colors =ann_colors)
dev.off()





#### Diffusion v1 #####
# Read-in object
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Diffusion pseudotime calculation
dpt <- DPT(dm)
es$pseudotime_dpt <- rank(dpt$dpt)
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

saveRDS(dm,paste0(path.ct,'data_',filename,'_dm.rds'))

# Plots
pdf(paste0(path.ct,'plots_',filename,'_diffusion.pdf'))

# UMAPs
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T,cols = colSTKvsNonSTKmain)
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T,cols = colSTKvsNonSTKmain) + NoLegend()
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T,cols = colPrimvsBMmain)
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T,cols = colPrimvsBMmain) + NoLegend()
DimPlot(seu, group.by = 'orig.ident', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'orig.ident', raster = T, shuffle = T) + NoLegend()

# 3D diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSTKvsNonSTKmain)
tmp_var<-seu$STK11mut_vs_NonSTK11mut
plot(dm, col = as.factor(tmp_var), main = 'STK11mut_vs_NonSTK11mut', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colPrimvsBMmain2)
tmp_var<-seu$PRIMARY_vs_BRAIN_METS
plot(dm, col = as.factor(tmp_var), main = 'PRIMARY_vs_BRAIN_METS', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette('ggplot2')
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette('ggplot2')
tmp_var<-seu$orig.ident
plot(dm, col = as.factor(tmp_var), main = 'orig.ident', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<- seu$pseudotime_dpt
plot(dm, col = sig_col, main = 'pseudotime_dpt', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)),
            posx = c(0.88, 0.9), posy = c(0, 0.65))

sig_col<-seu@assays$RNA@data['COL15A1',]
plot(dm, col = sig_col, main = 'COL15A1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

# 2D diffusion maps
set.seed(1)
dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))), 1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('cell_type_fine') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'STK11mut_vs_NonSTK11mut'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('STK11mut_vs_NonSTK11mut') + 
  theme(legend.title = element_blank())+
  scale_color_manual(values = colSTKvsNonSTKmain)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'PRIMARY_vs_BRAIN_METS'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('PRIMARY_vs_BRAIN_METS') + 
  theme(legend.title = element_blank())+
  scale_color_manual(values = colPrimvsBMmain2)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'orig.ident'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('orig.ident') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'pseudotime_dpt'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('pseudotime_dpt') + 
  theme(legend.title = element_blank()) + scale_color_viridis()

FeaturePlot(seu,features = paste0('DC',seq(1,9)),order = T,raster = F,pt.size = 0.1)

dev.off()


#### Violin plots #####
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
bm_comp<-subset(seu, orig.ident %notin% c('N561'))


genes <- c('VEGFC', 'DLL4', 'EFNB2', #arterial
           'RGCC','KIAA1217','ARHGAP18', #capillary
           'MMRN1', 'FLT4','SEMA3D','RELN', #lymphatic
           'IL1R1','NR2F2','CDH11', #venous
           'MKI67'
           )

pdf(paste0(path.ct,'plots_',filename,'_violin.pdf'))
VlnPlot(seu, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'none')+
  ggtitle('EC cell_type_fine (all samples)')

VlnPlot(seu, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T,split.by = 'STK11mut_vs_NonSTK11mut',
        cols = colSTKvsNonSTKmain) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom')+
  ggtitle('EC STK11mut_vs_NonSTK11mut (all samples)')

VlnPlot(bm_comp, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T,split.by = 'PRIMARY_vs_BRAIN_METS',
        cols = colPrimvsBMmain) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom')+
  ggtitle('EC PRIMARY_vs_BRAIN_METS (N561 excluded)')
dev.off()


#### Cell type frequencies ####
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
bm_comp<-subset(seu, orig.ident %notin% c('N561'))

df <- data.frame(orig.ident = seu$orig.ident, 
                 STK11mut_vs_NonSTK11mut = seu$STK11mut_vs_NonSTK11mut, 
                 cell_type_fine = seu$cell_type_fine, 
                 PRIMARY_vs_BRAIN_METS = seu$PRIMARY_vs_BRAIN_METS)

df_bm <- data.frame(orig.ident = bm_comp$orig.ident, 
                 STK11mut_vs_NonSTK11mut = bm_comp$STK11mut_vs_NonSTK11mut, 
                 cell_type_fine = bm_comp$cell_type_fine, 
                 PRIMARY_vs_BRAIN_METS = bm_comp$PRIMARY_vs_BRAIN_METS)

# cell_type_fine
df_fine = df %>%
  group_by(orig.ident, cell_type_fine, STK11mut_vs_NonSTK11mut) %>%
  tally() %>%
  group_by(orig.ident) %>%
  mutate(freq = n/sum(n))

df_fine_bm = df_bm %>%
  group_by(orig.ident, cell_type_fine, PRIMARY_vs_BRAIN_METS) %>%
  tally() %>%
  group_by(orig.ident) %>%
  mutate(freq = n/sum(n))


## Generate boxplots of cell-type frequencies
pdf(file = paste0(path.ct,'plots_',filename,'_celltype_freq.pdf'))
# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'STK11mut_vs_NonSTK11mut', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
  stat_compare_means(aes(group = STK11mut_vs_NonSTK11mut),label = 'p.format',
                     method = 'wilcox.test', size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ggtitle('EC cell_type_fine (all samples)') +
  scale_color_manual(values = colSTKvsNonSTKmain)

ggboxplot(df_fine_bm, x = 'cell_type_fine', y = 'freq', color = 'PRIMARY_vs_BRAIN_METS', add = 'jitter', 
          order = sort(unique(df_fine_bm$cell_type_fine))) + ylim(0, 1) + 
  stat_compare_means(aes(group = PRIMARY_vs_BRAIN_METS),label = 'p.format',
                     method = 'wilcox.test', size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ggtitle('EC cell_type_fine (N561 excluded)') +
  scale_color_manual(values = colPrimvsBMmain)

ggplot(seu@meta.data, aes(x = orig.ident, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('EC cell_type_fine') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x = STK11mut_vs_NonSTK11mut, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('EC cell_type_fine (all samples)') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(bm_comp@meta.data, aes(x = PRIMARY_vs_BRAIN_METS, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('EC cell_type_fine (N561 excluded)') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(0, 8, 0, 0, 'cm'))

dev.off()



