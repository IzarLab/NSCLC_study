# Author: Somnath Tagore, Ph.D. Title: T-cells volcano
# Script Name: volcano_tcells.R 
# Last Updated: 05/20/2021

# Packages required for this analysis
library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
library(hypeR)
library(DropletUtils)
'%notin%' <- Negate('%in%')

colExp<-c('#C82F71','#7473A6')


### Compare all expanded CD8+ T cells vs non-expanded CD8+ T cells
type <- 'Hyperexpanded vs Large/Medium'

# Read-in integrated object
patient <- readRDS('NSCLC_44_samples_T_cells.rds')

seu <- readRDS('NSCLC_44_samples_T_cells_v3.rds')
seu <- subset(seu, #sequencing == 'Single cells' &
                manual.annot == 'CD8+' &
                cloneType %in% c('Hyperexpanded (0.1 < X <= 1)', 'Large (0.01 < X <= 0.1)', 'Medium (0.001 < X <= 0.01)'))# &
                #both_chains == 'both')
seu <- ScaleData(seu)


#### Remove TRAV/TRBV genes ####
DefaultAssay(seu) <- 'RNA'
seu <- seu[grep('^TRAV|^TRBV', rownames(seu), value = T, invert = T), ]
seu <- NormalizeData(seu)


#### DEG ####
Idents(seu) <- seu$clone_size
markers <- FindMarkers(seu, ident.1 = 'Hyperexpanded (0.1 < X <= 1)', ident.2 = c('Large (0.01 < X <= 0.1)', 'Medium (0.001 < X <= 0.01)'),
                       min.pct = 0, logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA')
markers$gene <- rownames(markers)

tops <- markers %>%
  filter(p_val_adj < 0.05) %>%
  top_n(n = 10, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>%
                           filter(p_val_adj < 0.05) %>%
                           top_n(n = -10, wt = avg_log2FC))
markers$diffexpressed <- ifelse(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, 'Up',
                                ifelse(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05,
                                       'Down', 'NS'))
markers$label <- ifelse(markers$gene %in% tops$gene | markers$gene %in% c('TOX', 'TCF7'),
                        markers$gene, NA)

#ifelse(!dir.exists(file.path(paste0('data/DEG/clone_size/', type))),
#       dir.create(file.path(paste0('data/DEG/clone_size/', type)), recursive = T), FALSE)
write.csv(markers, paste0('type_markers_NSCLC_sc_', type, '_clone_size.csv'), row.names = F)


# Plot
pdf(paste0('type_plots_volcano_NSCLC_sc_', type, '_clone_size.pdf'))
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) +
  geom_point() +
  geom_text_repel() +
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) +
  labs(color = 'Log2FC\n(Hyperexpanded vs Large/Medium\n)') +
  theme_bw() +
  coord_cartesian(xlim = c(-3, 3)) +
  ggtitle(paste0('Hyperexpanded vs Large/Medium ', type)) +
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))

ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) +
  geom_point() +
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) +
  labs(color = 'Log2FC\n(Hyperexpanded vs Large/Medium\n)') +
  theme_bw() +
  coord_cartesian(xlim = c(-3, 3)) +
  ggtitle(paste0('Hyperexpanded vs Large/Medium ', type)) +
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
dev.off()
