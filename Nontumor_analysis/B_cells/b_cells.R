#!/usr/bin/env Rscript

### title: Analysis of B cells from NSCLC
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

path.ct<-'data/b_cells/'
filename<-'NSCLC_b_cells'

setwd('~/Documents/NSCLC/')

system(paste0("aws s3 sync s3://nscl-seq/Seurat/integrated_data/44_samples/nontumor_jana/b_cells/ ~/Documents/NSCLC/data/b_cells/ --exclude '*' --include 'Nontumor_group2_after_manual_annotation_v5*' --quiet"))

seu<- readRDS('data/b_cells/Nontumor_group2_after_manual_annotation_v5.rds')

# get module scores for sigs
sigs<-read.csv('~/brain_mets/signatures/bcell_signatures.csv',na.strings = '')
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
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.1 %in% c(2), 'Activated B cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.1 %in% c(3), 'Cycling cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.1 %in% c(0,1,4,5,7), 'Plasma cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.1 %in% c(6), 'Stromal cells', seu$cell_type_fine)

seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Stromal cells'), 'Stromal cells', 'B/Plasma cells')


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
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'patient',shuffle = T,raster = T)+NoLegend())
DimPlot(seu,group.by = 'cell_type_fine',shuffle = T,raster = T)
FeatureScatter(seu,'G2M.Score','S.Score',group.by = 'cell_cycle', pt.size = 0.01)
print(DimPlot(seu, reduction = "umap",label = F,group.by = 'cell_cycle',shuffle = T,raster = T))

DimPlot(seu, reduction = "umap",label = T,group.by = 'manual.annot',repel = T,
        label.size = 3,shuffle = T,raster = T)

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_bped_fine',repel = T,
        label.size = 3,shuffle = T,raster = T)+NoLegend()

DimPlot(seu, reduction = "umap",label = T,group.by = 'celltype_hpca_main',repel = T,
        label.size = 2.5,shuffle = T,raster = T)

DotPlot(seu,assay = 'RNA',features = c('CD27','SLAMF7','IRF4','SDC1','CD38','PRDM1','XBP1', #plasma
                                       'CD24','MS4A1','CD19', #b cell
                                       'CD22','IGHD','CCR6','CR2', #memory
                                       'TCL1A','LAIR1','IGHM', #naive
                                       'IL2RA','BLK', #activated
                                       'MKI67'
),
group.by = 'integrated_snn_res.0.1', dot.scale = 5)+ scale_color_viridis() + coord_flip()+
  theme(axis.text.y = element_text(size=7))

DotPlot(seu,assay = 'RNA',features = c('CD27','SLAMF7','IRF4','SDC1','CD38','PRDM1','XBP1', #plasma
                                       'CD24','MS4A1','CD19', #b cell
                                       'CD22','IGHD','CCR6', #memory
                                       'IL2RA','BLK', #activated
                                       'TCL1A','LAIR1','IGHM','IGHD', #naive
                                       'MKI67'
),
group.by = 'cell_type_fine', dot.scale = 5)+ scale_color_viridis() + coord_flip()+ RotatedAxis()+
  theme(axis.text.y = element_text(size=7))

DotPlot(seu,assay = 'RNA',features = paste0(names(sigs),1),
        group.by = 'integrated_snn_res.0.1', dot.scale = 3)+ 
  scale_color_viridis() + coord_flip()


# naive
print(FeaturePlot(seu, features = c('rna_TCL1A','rna_LAIR1','rna_IGHM'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))

# memory
print(FeaturePlot(seu, features = c('rna_CD22','rna_IGHD','rna_CCR6','rna_CR2'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))

# activated
print(FeaturePlot(seu, features = c('rna_IL2RA','rna_BLK','rna_CD24','rna_MS4A1'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))

# plasma
print(FeaturePlot(seu, features = c('rna_JCHAIN','rna_MZB1','rna_XBP1','rna_CD38'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = T))


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
ct <- seu@meta.data %>% select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct,'data_',filename,'_celltype.csv'))


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


#### Diffusion #####
# Read-in object and subset to B cells
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
#seu <- subset(seu, cell_type_main %in% c('B/Plasma cells'))
seu <- subset(seu, cell_type_fine %in% c('Activated B cells', 'Plasma cells'))
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
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T)
DimPlot(seu, group.by = 'STK11mut_vs_NonSTK11mut', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T)
DimPlot(seu, group.by = 'PRIMARY_vs_BRAIN_METS', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'orig.ident', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'orig.ident', raster = T, shuffle = T) + NoLegend()

# 3D diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette("ggplot2")
tmp_var<-seu$STK11mut_vs_NonSTK11mut
plot(dm, col = as.factor(tmp_var), main = 'STK11mut_vs_NonSTK11mut', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

#palette(colBP)
tmp_var<-seu$PRIMARY_vs_BRAIN_METS
plot(dm, col = as.factor(tmp_var), main = 'PRIMARY_vs_BRAIN_METS', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

#palette(colBcell)
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

#palette(colBcell)
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

sig_col<-seu@assays$RNA@data['XBP1',]
plot(dm, col = sig_col, main = 'XBP1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

# 2D diffusion maps
set.seed(1)
dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))), 1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('cell_type_fine') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'STK11mut_vs_NonSTK11mut'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('STK11mut_vs_NonSTK11mut') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'PRIMARY_vs_BRAIN_METS'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('PRIMARY_vs_BRAIN_METS') + 
  theme(legend.title = element_blank())

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
seu <- subset(seu, cell_type_main %in% c('B/Plasma cells'))
bm_comp<-subset(seu, orig.ident %notin% c('N561'))


genes <- c('BANK1', 'BLK', # activated B
           'MKI67','TOP2A',
           'SDC1', 'CD38', 'PRDM1' # plasma
           )

pdf('data/b_cells/plots_NSCLC_b_cells_violin.pdf')
VlnPlot(seu, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'none')+
  ggtitle('B cells cell_type_fine (all samples)')

VlnPlot(seu, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T,split.by = 'STK11mut_vs_NonSTK11mut') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom')+
  ggtitle('B cells STK11mut_vs_NonSTK11mut (all samples)')

VlnPlot(bm_comp, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T,split.by = 'PRIMARY_vs_BRAIN_METS') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom')+
  ggtitle('B cells PRIMARY_vs_BRAIN_METS (N561 excluded)')
dev.off()


#### Cell type frequencies ####
seu <- readRDS(paste0(path.ct,'data_',filename,'.rds'))
seu <- subset(seu, cell_type_main %in% c('B/Plasma cells'))
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
pdf(file = 'data/b_cells/plots_NSCLC_celltype_freq.pdf')
# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'STK11mut_vs_NonSTK11mut', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
  stat_compare_means(aes(group = STK11mut_vs_NonSTK11mut),label = 'p.format',
                     method = 'wilcox.test', size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ggtitle('B cells cell_type_fine (all samples)') 

ggboxplot(df_fine_bm, x = 'cell_type_fine', y = 'freq', color = 'PRIMARY_vs_BRAIN_METS', add = 'jitter', 
          order = sort(unique(df_fine_bm$cell_type_fine))) + ylim(0, 1) + 
  stat_compare_means(aes(group = PRIMARY_vs_BRAIN_METS),label = 'p.format',
                     method = 'wilcox.test', size = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ggtitle('B cells cell_type_fine (N561 excluded)') 

ggplot(seu@meta.data, aes(x = orig.ident, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('B cells cell_type_fine') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(seu@meta.data, aes(x = STK11mut_vs_NonSTK11mut, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('B cells cell_type_fine (all samples)') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(bm_comp@meta.data, aes(x = PRIMARY_vs_BRAIN_METS, fill = cell_type_fine)) + 
  geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
  ylab('Cell type fraction') + ggtitle('B cells cell_type_fine (N561 excluded)') + 
  theme(legend.position = 'right', axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        plot.margin = margin(0, 8, 0, 0, 'cm'))

dev.off()



