# Author: Somnath Tagore, Ph.D. 
# Title: Cluster 21 analysis 
# Script Name: cluster_21_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs analysis on Cluster 21
#
# -------------------------
#!/usr/bin/env Rscript

library(matrixStats)
library(copynumber)
library(infercnv)
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(gplots)
library(singscore)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(reshape2)
library(ggpubr)

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

# 1) Cluster 21 identification
# -----------------------------

pat <- 'Integrated'

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

# read the sample-specific tumor data
PA141 <- readRDS(file="PA141_seurat.rds")
PA141
PA141.tumor.labels <- read.table(file="PA141_tumor_most_upd.csv")
row.names(PA141.tumor.labels) <- PA141.tumor.labels[,1]
PA141.tumor.labels <- PA141.tumor.labels[,-1]
PA141[['tumor_most_upd']]<-ifelse(colnames(PA141) %in% rownames(PA141.tumor.labels),'Tumor','Non-Tumor')
PA141.tumor<-subset(PA141, tumor_most_upd == 'Tumor')
PA141.tumor@meta.data$orig.ident <- rep('PA141',length(colnames(PA141.tumor)))

# read all other 43 samples

# merge samples
nsclc.combined <- merge(PA141.tumor, y =c(PA125.tumor,PA104.tumor,PA080.tumor,PA076.tumor,PA072.tumor,PA070.tumor,PA068.tumor,PA067.tumor,PA060.tumor,PA056.tumor,PA054.tumor,PA048.tumor,PA043.tumor,PA042.tumor,PA034.tumor,PA025.tumor,PA019.tumor,PA005.tumor,PA004.tumor,PA001.tumor,N586.tumor,N561.tumor,N254.tumor,STK_5dot2.tumor,STK_5dot1.tumor,STK_3.tumor,STK_2.tumor,STK_22dot2.tumor,STK_21.tumor,STK_20.tumor,STK_1.tumor,STK_18.tumor,STK_15.tumor,STK_14.tumor,KRAS_8.tumor,KRAS_7.tumor,KRAS_6.tumor,KRAS_4.tumor,KRAS_17.tumor,KRAS_13.tumor,KRAS_12.tumor,KRAS_11.tumor,KRAS_10.tumor), add.cell.ids = c("PA141", "PA125","PA104","PA080","PA076","PA072","PA070","PA068","PA067","PA060","PA056","PA054","PA048","PA043","PA042","PA034","PA025","PA019","PA005","PA004","PA001","N586","N561","N254","STK_5dot2","STK_5dot1","STK_3","STK_2","STK_22dot2","STK_21","STK_20","STK_1","STK_18","STK_15","STK_14","KRAS_8","KRAS_7","KRAS_6","KRAS_4","KRAS_17","RAS_13","KRAS_12","KRAS_11","KRAS_10"), project = "NSCLC")
nsclc.combined

table(nsclc.combined$orig.ident)

immune.combined<-nsclc.combined
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "RNA"

# Run the standard workflow for visualization and clustering
immune.combined <- NormalizeData(object = immune.combined)
immune.combined <- FindVariableFeatures(immune.combined, do.plot = F, display.progress = F)
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined,npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

pdf(file = "NSCLC_KRAS_STK_samples_tumor_updated_overlay_44_samples_merge_no_anchor_test.pdf",height = 20,width=20)

# Visualization
print(DimPlot(immune.combined, reduction = "pca"))
print(DimPlot(immune.combined, reduction = "umap"))
#print(DimPlot(immune.combined, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(immune.combined, reduction = "pca", group.by = "orig.ident", cols = c('PA001' = '#808080', 'PA004' = '#d3d3d3', 'PA005' = '#2f4f4f',
                                                                                    'PA019' = '#556b2f', 'PA025' = '#8b4513','PA034' = '#2e8b57', 'PA042' = '#228b22',
                                                                                    'PA043' = '#7f0000','PA048' = '#191970', 'PA054' = '#808000', 'PA056' = '#b8860b',
                                                                                    'PA060' = '#008b8b',
                                                                                    'PA067' = '#4682b4', 'PA068' = '#d2691e','PA070' = '#9acd32','PA072' = '#cd5c5c',
                                                                                    'PA076' = '#00008b', 'PA080' = '#32cd32','PA104' = '#8fbc8f','PA125' = '#8b008b',
                                                                                    'PA141' = '#b03060','N254' = '#ff4500','N561' = '#00ced1','N586' = '#ffa500',
                                                                                    'KRAS_10' = '#ffd700', 'KRAS_11' = '#6a5acd', 'KRAS_12' = '#deb887',
                                                                                    'KRAS_13' = '#00ff00', 'KRAS_17' = '#00fa9a', 'KRAS_4' = '#dc143c',
                                                                                    'KRAS_6' = '#0000ff', 'KRAS_7' = '#a020f0', 'KRAS_8' = '#adff2f',
                                                                                    'STK_1' = '#da70d6', 'STK_14' = '#ff00ff', 'STK_15' = '#1e90ff',
                                                                                    'STK_18' = '#f0e68c', 'STK_2' = '#dda0dd','STK_20' = '#90ee90',
                                                                                    'STK_21' = '#ffa07a', 'STK_22dot2' = '#87cefa', 'STK_3' = '#7fffd4',
                                                                                    'STK_5dot1' = '#ff69b4','STK_5dot1' = '#ffb6c1'),label = TRUE, repel = TRUE)+ labs(title = paste0("PCA All 44 Samples NSCLC Tumor Cells (Merge, no integration)")))


#print(DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE))
#print(DimPlot(immune.combined, reduction = "umap", group.by = "ident", label = TRUE, repel = TRUE))
print(DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident", cols = c('PA001' = '#808080', 'PA004' = '#d3d3d3', 'PA005' = '#2f4f4f',
                                                                                     'PA019' = '#556b2f', 'PA025' = '#8b4513','PA034' = '#2e8b57', 'PA042' = '#228b22',
                                                                                     'PA043' = '#7f0000','PA048' = '#191970', 'PA054' = '#808000', 'PA056' = '#b8860b',
                                                                                     'PA060' = '#008b8b',
                                                                                     'PA067' = '#4682b4', 'PA068' = '#d2691e','PA070' = '#9acd32','PA072' = '#cd5c5c',
                                                                                     'PA076' = '#00008b', 'PA080' = '#32cd32','PA104' = '#8fbc8f','PA125' = '#8b008b',
                                                                                     'PA141' = '#b03060','N254' = '#ff4500','N561' = '#00ced1','N586' = '#ffa500',
                                                                                     'KRAS_10' = '#ffd700', 'KRAS_11' = '#6a5acd', 'KRAS_12' = '#deb887',
                                                                                     'KRAS_13' = '#00ff00', 'KRAS_17' = '#00fa9a', 'KRAS_4' = '#dc143c',
                                                                                     'KRAS_6' = '#0000ff', 'KRAS_7' = '#a020f0', 'KRAS_8' = '#adff2f',
                                                                                     'STK_1' = '#da70d6', 'STK_14' = '#ff00ff', 'STK_15' = '#1e90ff',
                                                                                     'STK_18' = '#f0e68c', 'STK_2' = '#dda0dd','STK_20' = '#90ee90',
                                                                                     'STK_21' = '#ffa07a', 'STK_22dot2' = '#87cefa', 'STK_3' = '#7fffd4',
                                                                                     'STK_5dot1' = '#ff69b4','STK_5dot1' = '#ffb6c1'),label = TRUE, repel = TRUE)+ labs(title = paste0("UMAP All 44 Samples NSCLC Tumor Cells (Merge, no integration)")))


saveRDS(immune.combined, file = 'NSCLC_tumor_just_merge_no_anchor.rds')

# 2) Cluster 21 Quality Control
# -----------------------------
#read metadata
ichorcna_data <- read.csv(file="NSCLC_tumor_just_merge_no_anchor_metadata.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Cluster_21_vs_rest,
                              percent.mt,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$Cluster_21_vs_rest <- factor(data_ordered$Cluster_21_vs_rest,
                                          levels = levels(group_ordered))
# MT %age: Cluster 21 vs Rest
yplot <- ggboxplot(data_ordered, x = "Cluster_21_vs_rest", y = "percent.mt", title = "MT %age: Cluster 21 vs Rest",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Cluster_21_vs_rest", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 42))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Cluster_21_vs_rest)),
           y = aggregate(percent.mt ~ Cluster_21_vs_rest, data_ordered, median)[ , 2],
           label = table(data_ordered$Cluster_21_vs_rest),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00')) + scale_fill_manual(values=c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00'))
ggsave("Cluster_21_vs_Rest_MT.pdf",height = 5,width = 5)
#dev.off()

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Cluster_21_vs_rest,
                              nCount_RNA,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$Cluster_21_vs_rest <- factor(data_ordered$Cluster_21_vs_rest,
                                          levels = levels(group_ordered))

# UMI counts: Cluster 21 vs Rest
yplot <- ggboxplot(data_ordered, x = "Cluster_21_vs_rest", y = "nCount_RNA", title = "UMI counts: Cluster 21 vs Rest",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Cluster_21_vs_rest", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 42))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Cluster_21_vs_rest)),
           y = aggregate(nCount_RNA ~ Cluster_21_vs_rest, data_ordered, median)[ , 2],
           label = table(data_ordered$Cluster_21_vs_rest),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00')) + scale_fill_manual(values=c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00'))
ggsave("Cluster_21_vs_Rest_UMI.pdf",height = 5,width = 5)
#dev.off()

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Cluster_21_vs_rest,
                              nFeature_RNA,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$Cluster_21_vs_rest <- factor(data_ordered$Cluster_21_vs_rest,
                                          levels = levels(group_ordered))

# Gene counts: Cluster 21 vs Rest
yplot <- ggboxplot(data_ordered, x = "Cluster_21_vs_rest", y = "nFeature_RNA", title = "Gene counts: Cluster 21 vs Rest",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Cluster_21_vs_rest", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 42))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Cluster_21_vs_rest)),
           y = aggregate(nFeature_RNA ~ Cluster_21_vs_rest, data_ordered, median)[ , 2],
           label = table(data_ordered$Cluster_21_vs_rest),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00')) + scale_fill_manual(values=c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00'))
ggsave("Cluster_21_vs_Rest_Gene_Count.pdf",height = 5,width = 5)
#dev.off()

# 3) CIN70 signature Cluster_21 Prim vs BM
# ----------------------------------------
#read metadata
ichorcna_data <- read.csv(file="Cluster_21_NSCLC_meta_data.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
table(ichorcna_data$PRIMARY_vs_BRAIN_METS)
table(ichorcna_data$STK11_mut_vs_Non_STK11_mut)
ichorcna_data[1:10,]

yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CIN70_signature_score_upd", title = "Cluster_21 - CIN70 Signature: PRIMARY_vs_BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())

yplot + stat_compare_means()

ggsave("Cluster_21_CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)
dev.off()

# 4) CIN70 signature Cluster_21 vs rest
# -------------------------------------

#read metadata
ichorcna_data <- read.csv(file="CIN70_cluster21_vs_rest.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Cluster_21_vs_rest,
                              Cluster_21_vs_rest_CIN70,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$Cluster_21_vs_rest <- factor(data_ordered$Cluster_21_vs_rest,
                                          levels = levels(group_ordered))

yplot <- ggboxplot(data_ordered, x = "Cluster_21_vs_rest", y = "Cluster_21_vs_rest_CIN70", title = "CIN70 Score: Cluster 21 vs Rest",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Cluster_21_vs_rest", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 42))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Cluster_21_vs_rest)),
           y = aggregate(Cluster_21_vs_rest_CIN70 ~ Cluster_21_vs_rest, data_ordered, median)[ , 2],
           label = table(data_ordered$Cluster_21_vs_rest),
           col = "red",
           vjust = - 1)+ stat_compare_means()#+ scale_color_manual(values = c('ALK_MUT'='#00CED1','EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082')) + scale_fill_manual(values=c('ALK_MUT'='#00CED1','EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082'))
ggsave("Cluster_21_vs_Rest_CIN70_signature.pdf",height = 5,width = 5)
#dev.off()

# 5) CIN70 signature Cluster_21 STK11_mut_vs_STK11_wt
# ---------------------------------------------------

#read metadata
ichorcna_data <- read.csv(file="Cluster_21_NSCLC_meta_data.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "CIN70_signature_score_upd", title = "Cluster_21 - CIN70 Signature: STK11_mut_vs_Non_STK11_mut",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
#yplot + stat_compare_means(comparisons = my_comparisons)
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/Ben_Izar_Project/CIN_signature/Cluster_21_CIN70_signature_STK11_mut_vs_Non_STK11_mut.pdf",height = 10,width = 10)
dev.off()

# 6) Cluster_21_ signature
# ------------------------

cols = c(
  'PA001' = '#808080', 'PA004' = '#d3d3d3', 'PA005' = '#2f4f4f',
  'PA019' = '#556b2f', 'PA025' = '#8b4513','PA034' = '#2e8b57', 'PA042' = '#228b22',
  'PA043' = '#7f0000','PA048' = '#191970', 'PA054' = '#808000', 'PA056' = '#b8860b',
  'PA060' = '#008b8b',
  'PA067' = '#4682b4', 'PA068' = '#d2691e','PA070' = '#9acd32','PA072' = '#cd5c5c',
  'PA076' = '#00008b', 'PA080' = '#32cd32','PA104' = '#8fbc8f','PA125' = '#8b008b',
  'PA141' = '#b03060','N254' = '#ff4500',#'N561' = '#00ced1',
  'N586' = '#ffa500',
  'KRAS_10' = '#ffd700', 'KRAS_11' = '#6a5acd', 'KRAS_12' = '#deb887',
  'KRAS_13' = '#00ff00', 'KRAS_17' = '#00fa9a', 'KRAS_4' = '#dc143c',
  'KRAS_6' = '#0000ff', 'KRAS_7' = '#a020f0', 'KRAS_8' = '#adff2f',
  'STK_1' = '#da70d6', 'STK_14' = '#ff00ff', 'STK_15' = '#1e90ff',
  'STK_18' = '#f0e68c', 'STK_2' = '#dda0dd','STK_20' = '#90ee90',
  'STK_21' = '#ffa07a', 'STK_22dot2' = '#87cefa', 'STK_3' = '#7fffd4',
  'STK_5dot1' = '#ff69b4','STK_5dot1' = '#ffb6c1'
)

group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(orig.ident,
                              desc(cluster_21_signature1),
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$orig.ident <- factor(data_ordered$orig.ident,
                                  levels = levels(group_ordered))

# Cluster_21 signature: All samples outside cluster_21
yplot <- ggboxplot(data_ordered, x = "orig.ident", y = "cluster_21_signature1", title = "Cluster_21 signature: All samples outside cluster_21",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "orig.ident", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 24))+theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$orig.ident)),
           y = aggregate(cluster_21_signature1 ~ orig.ident, data_ordered, median)[ , 2],
           label = table(data_ordered$orig.ident),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = cols) + scale_fill_manual(values=cols)

ggsave("Cluster_21_signature_all_samples_outside_cluster_21_v2.pdf",height = 15,width = 25)

# 7) Cluster_21 signature: Non-tumor groups
# --------------------------------------
ichorcna_data <- read.csv(file="NSCLC_tumor_nontumor_just_merge_no_anchor_cluster_21_signature_updated_1.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(Non_tumor_groups,
                              desc(cluster21_signature),
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$Non_tumor_groups <- factor(data_ordered$Non_tumor_groups,
                                        levels = levels(group_ordered))

yplot <- ggboxplot(data_ordered, x = "Non_tumor_groups", y = "cluster21_signature", title = "Cluster_21 signature: Non-tumor groups",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Non_tumor_groups", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 24))+theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Non_tumor_groups)),
           y = aggregate(cluster21_signature ~ Non_tumor_groups, data_ordered, median)[ , 2],
           label = table(data_ordered$Non_tumor_groups),
           col = "red",
           vjust = - 1)+ stat_compare_means()#+ scale_color_manual(values = cols) + scale_fill_manual(values=cols)

ggsave("Cluster_21_signature_non_tumor_v2.pdf",height = 7,width = 7)

# 8) LISI Score: Cluster_21
# -------------------------

# LISI calculation ???

ichorcna_data_bi <- read.csv(file="~/Documents/Ben_Izar_Project/cluster_21/lisi.cluster_21_others.csv")
#ichorcna_data <- as.data.frame(ichorcna_data)
ichorcna_data_bi <- as.data.frame(ichorcna_data_bi)
table(ichorcna_data_bi$Cluster_21_vs_Others)
yplot <- ggboxplot(ichorcna_data_bi, x = "Cluster_21_vs_Others", y = "mycluster", title = "Cluster_21 vs Not_Cluster_21: LISI Score",
                   color = "Cluster_21_vs_Others", fill = "Cluster_21_vs_Others", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 16))+   # Add counts by group to boxplot
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(ichorcna_data_bi$Cluster_21_vs_Others)),
           y = aggregate(mycluster ~ Cluster_21_vs_Others, ichorcna_data_bi, median)[ , 2],
           label = table(ichorcna_data_bi$Cluster_21_vs_Others),
           col = "red",
           vjust = - 1)
ggsave("lisi.cluster_21_others.pdf",height = 7,width = 7)

# 9) Cluster_21 stress signature
# ------------------------------

#stress signatures
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/cluster_21_stress_denisenko.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(Cluster_21_vs_rest,
                              stress_denisenko,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$Cluster_21_vs_rest <- factor(data_ordered$Cluster_21_vs_rest,
                                          levels = levels(group_ordered))

# Stress module: Cluster 21 vs Rest
# Stress brink: Cluster 21 vs Rest
# Stress vanhove: Cluster 21 vs Rest
# Stress marsh: Cluster 21 vs Rest
# Stress denisenko: Cluster 21 vs Rest 

yplot <- ggboxplot(data_ordered, x = "Cluster_21_vs_rest", y = "stress_denisenko", title = "Stress denisenko: Cluster 21 vs Rest",
                   # color = "STK11_MUT_vs_Other_genes", 
                   fill = "Cluster_21_vs_rest", #palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"))#+ scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                         values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) 
#yplot + stat_compare_means()+theme(text = element_text(size = 20))
yplot + stat_compare_means()+theme(text = element_text(size = 42))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Cluster_21_vs_rest)),
           y = aggregate(stress_denisenko ~ Cluster_21_vs_rest, data_ordered, median)[ , 2],
           label = table(data_ordered$Cluster_21_vs_rest),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00')) + scale_fill_manual(values=c('Cluster_21'='#c00000','Not_Cluster_21'='#CCFF00'))
ggsave("Cluster_21_vs_Rest_stress_denisenko.pdf",height = 5,width = 5)

# 10) cluster_21 top300 genes signature on nontumor samples
# --------------------------------------------------------

patients.integrated <- readRDS(file = "Integrated_NSCLC_KRAS_STK_samples_nontumor_44_samples_after_manual_annotation_update_1.rds")

cluster21_top300_genes <- read.csv(file="cluster21_top300_genes.csv")

manual.annot.refined.nontumors.labels <- read.csv(file="manual.annot.refined.nontumors.csv")
patients.integrated[['manual.annot.refined.nontumors']]<-manual.annot.refined.nontumors.labels[,1]

#DefaultAssay(patients.integrated) <- "RNA"
DefaultAssay(patients.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
#patients.integrated <- NormalizeData(patients.integrated)
#patients.integrated <- FindVariableFeatures(patients.integrated)
patients.integrated <- ScaleData(patients.integrated, verbose = FALSE)
patients.integrated <- RunPCA(patients.integrated, npcs = 30, verbose = FALSE)
patients.integrated <- RunUMAP(patients.integrated, reduction = "pca", dims = 1:30)
patients.integrated <- FindNeighbors(patients.integrated, reduction = "pca", dims = 1:30)

#DefaultAssay(patients.integrated) <- "RNA"
patients.integrated <- AddModuleScore(
  object = patients.integrated,
  features = list(cluster21_top300_genes[,1]),
  assay = 'RNA',
  ctrl = 5,
  name = 'cluster21_top300_genes'
)

patients.integrated[['cluster21_top300_genes']] <- CreateAssayObject(data = t(x = FetchData(object = patients.integrated, vars = 'cluster21_top300_genes1')))

whichsignature <- 'cluster21_top300_genes'

pdf(file = paste0("NSCLC_nontumor_cluster_21_signature_", whichsignature ,"_updated.pdf"),height = 20,width=20)

VlnPlot(patients.integrated,features = "cluster21_top300_genes1", pt.size = 0,assay = 'RNA',flip = T,split.by = 'manual.annot',split.plot = F)

dev.off()

# 11) Cluster 21 signature on CCLE
# --------------------------------
# '%notin%' <- Negate('%in%')
# 

#### DESeq2

gex.11 <- read.csv(file="CCLE_Lung_Cancer_Expression_22Q2_Public_subsetted.csv")
row.names(gex.11) <- gex.11[,1]
gex.11<-gex.11[,-1]
gex.11 <- t(gex.11)

dim(gex.11)

gex.11[1:5,1:5]

cts <- as.matrix(gex.11)

# Apply SingScore on bulk data
rankData <- rankGenes(cts)
dim(rankData)

#colnames(rankData) = make.names(colnames(rankData), unique=TRUE)

#write.csv(t(rankData),file="rankData.LUAD.csv")
#write.csv(t(rankData),file="rankData.LUSC.csv")
rankData[1:5,1:5]
#markers_cin <- read.csv('CIN70_signature.csv')
markers_cluster_21 <- read.csv('~/Documents/Ben_Izar_Project/cluster_21/cluster_21_signature.csv')

scoredf <- simpleScore(rankData, upSet = markers_cluster_21[,1])
#write.csv(scoredf,file="singscore.LUAD.csv")
write.csv(scoredf,file="singscore.CCLE.lung.markers.cluster_21.csv")

# Cluster 21 signature CCLE
#TCGA LUAD All STK_mut vs All Non_STK_mut
library(ggplot2)
library(ggpubr)

ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/cluster_21/CCLE.lung.signature.cluster_21.updated.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

#CCLE (Lung) - Cluster_21_signature \n Primary vs Metastasi
yplot <- ggboxplot(ichorcna_data, x = "Primary_vs_Mets", y = "Cluster_21_signature", title = "CCLE (Lung) - Cluster_21_signature \n Primary vs Metastasis",
                   color = "Primary_vs_Mets", fill = "Primary_vs_Mets", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
#yplot  + stat_compare_means(method = "t.test")
ggsave("CCLE.lung.signature.cluster_21.Prim.vs.Mets.pdf",height = 15,width = 15)

#CCLE (Lung) - Cluster_21_signature \n Primary vs Metastasis
yplot <- ggboxplot(ichorcna_data, x = "CCLE_Cell_Line", y = "Cluster_21_signature", title = "CCLE (Lung) - Cluster_21_signature \n Primary vs Metastasis",
                   color = "CCLE_Cell_Line", fill = "CCLE_Cell_Line",# palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20),legend.position = "none")
ggsave("CCLE.lung.signature.cluster_21.CCLE_Cell_Line.pdf",height = 15,width = 15)

pdf("CCLE.lung.signature.cluster_21.CCLE_Cell_Line_1.pdf",height = 15,width = 15)
ggplot(data = ichorcna_data,
       mapping = aes(x = CCLE_Cell_Line, y = Cluster_21_signature)) +
       geom_bar(stat = "identity", position = "dodge") + 
       labs(x=NULL, y="Cluster_21_signature") +
       coord_flip() +
       theme_bw()+theme(text = element_text(size = 7)) + ggtitle("CCLE (Lung) - Cluster_21_signature \n All Lung Cancer Cell lines") 
dev.off()

# Cluster 21 signature: Public databases
# --------------------------------------

# PMID_32385277_cluster_21_signature
ichorcna_data <- read.csv(file="PMID_32385277_cluster_21_signature.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Tissue.origins,
                              # reorder(STK11_MUT_vs_Others,
                              Cluster_21_signature,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# data_ordered$STK11_MUT_vs_Others <- factor(data_ordered$STK11_MUT_vs_Others,
#                                            levels = levels(group_ordered))
# Create data with reordered group levels
data_ordered$Tissue.origins <- factor(data_ordered$Tissue.origins,
                                      levels = levels(group_ordered))

yplot <-   ggplot(data_ordered,                             # Draw grouped boxplot with default order
                  aes(#x = STK11_MUT_vs_Others,
                    x = Tissue.origins,
                    y = Cluster_21_signature,
                    fill = Tissue.origins)) + ggtitle(paste0("PMID: 32385277 - Cluster_21_signature"))+
  geom_boxplot() 

yplot + stat_compare_means()+theme(text = element_text(size = 24))+
  annotate("text",
           x = 1:length(table(data_ordered$Tissue.origins)),
           y = aggregate(Cluster_21_signature ~ Tissue.origins, data_ordered, median)[ , 2],
           label = table(data_ordered$Tissue.origins),
           col = "black",
           vjust = - 1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")
ggsave("PMID_32385277.cluster_21.signature.v2.pdf",height = 7,width = 7)

# PMID_34663877_cluster_21_signature
ichorcna_data <- read.csv(file="public_NSCLC_data/PMID_34663877_cluster_21_signature.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Tumor_stage,
                              # reorder(STK11_MUT_vs_Others,
                              Cluster_21_signature,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$Tumor_stage <- factor(data_ordered$Tumor_stage,
                                   levels = levels(group_ordered))

yplot <-   ggplot(data_ordered,                             # Draw grouped boxplot with default order
                  aes(#x = STK11_MUT_vs_Others,
                    x = Tumor_stage,
                    y = Cluster_21_signature,
                    fill = Tumor_stage)) + ggtitle(paste0("PMID: 34663877 - Cluster_21_signature"))+
  geom_boxplot() 

yplot + stat_compare_means()+theme(text = element_text(size = 24))+
  annotate("text",
           x = 1:length(table(data_ordered$Tumor_stage)),
           y = aggregate(Cluster_21_signature ~ Tumor_stage, data_ordered, median)[ , 2],
           label = table(data_ordered$Tumor_stage),
           col = "black",
           vjust = - 1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "none")+theme_bw()
ggsave("PMID_34663877.cluster_21.signature.v2.pdf",height = 10,width = 10)

# PMID_33953163_cluster_21_signature
ichorcna_data <- read.csv(file="PMID_33953163_cluster_21_signature.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(Lung_cancer,
                              # reorder(STK11_MUT_vs_Others,
                              Cluster_21_signature,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              
# Create data with reordered group levels

data_ordered$Lung_cancer <- factor(data_ordered$Lung_cancer,
                                   levels = levels(group_ordered))

yplot <-   ggplot(data_ordered,                             # Draw grouped boxplot with default order
                  aes(#x = STK11_MUT_vs_Others,
                    x = Lung_cancer,
                    y = Cluster_21_signature,
                    fill = Lung_cancer)) + ggtitle(paste0("PMID: 33953163 (Cluster_21_signature)"))+
  geom_boxplot() 

yplot + stat_compare_means()+theme(text = element_text(size = 24))+
  annotate("text",
           x = 1:length(table(data_ordered$Lung_cancer)),
           y = aggregate(Cluster_21_signature ~ Lung_cancer, data_ordered, median)[ , 2],
           label = table(data_ordered$Lung_cancer),
           col = "black",
           vjust = - 1)+theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("PMID_33953163.cluster_21.signature.v2.pdf",height = 10,width = 10)

# 12) Cluster 21: Sanity checks
# -----------------------------
pat <- 'Integrated'

NSCLC_tumor_just_merge_no_anchor <- readRDS(file="NSCLC_tumor_just_merge_no_anchor.rds")

nsclc.combined <- NSCLC_tumor_just_merge_no_anchor

table(nsclc.combined$orig.ident)

# Tumor overlay: Average CNV
infercnv.labels <- read.table(file=paste0("NSCLC_KRAS_STK_tumor_44_samples_infercnv_avg.csv"))

nsclc.combined[["avg_infercnv"]]<-infercnv.labels[,1]

for(i in 1:10)
{
  temp<-subset(nsclc.combined,subset=seurat_clusters==i)
  print(quantile(temp@meta.data$avg_infercnv))
}

for(i in 1:10)
{
  cl0<-subset(nsclc.combined,subset=seurat_clusters=='0')
  quantile(cl0@meta.data$avg_infercnv)
}
cl0<-subset(nsclc.combined,subset=seurat_clusters=='0')
quantile(cl0@meta.data$avg_infercnv)
cl1<-subset(nsclc.combined,subset=seurat_clusters=='1')
quantile(cl1@meta.data$avg_infercnv)
cl2<-subset(nsclc.combined,subset=seurat_clusters=='2')
quantile(cl2@meta.data$avg_infercnv)
cl3<-subset(nsclc.combined,subset=seurat_clusters=='3')
quantile(cl3@meta.data$avg_infercnv)
cl4<-subset(nsclc.combined,subset=seurat_clusters=='4')
quantile(cl4@meta.data$avg_infercnv)
cl5<-subset(nsclc.combined,subset=seurat_clusters=='5')
quantile(cl5@meta.data$avg_infercnv)
cl6<-subset(nsclc.combined,subset=seurat_clusters=='6')
quantile(cl6@meta.data$avg_infercnv)
cl7<-subset(nsclc.combined,subset=seurat_clusters=='7')
quantile(cl7@meta.data$avg_infercnv)
cl8<-subset(nsclc.combined,subset=seurat_clusters=='8')
quantile(cl8@meta.data$avg_infercnv)
cl9<-subset(nsclc.combined,subset=seurat_clusters=='9')
quantile(cl9@meta.data$avg_infercnv)
cl10<-subset(nsclc.combined,subset=seurat_clusters=='10')
quantile(cl10@meta.data$avg_infercnv)
cl11<-subset(nsclc.combined,subset=seurat_clusters=='11')
quantile(cl11@meta.data$avg_infercnv)
cl12<-subset(nsclc.combined,subset=seurat_clusters=='12')
quantile(cl12@meta.data$avg_infercnv)
cl13<-subset(nsclc.combined,subset=seurat_clusters=='13')
quantile(cl13@meta.data$avg_infercnv)
cl14<-subset(nsclc.combined,subset=seurat_clusters=='14')
quantile(cl14@meta.data$avg_infercnv)
cl15<-subset(nsclc.combined,subset=seurat_clusters=='15')
quantile(cl15@meta.data$avg_infercnv)
cl16<-subset(nsclc.combined,subset=seurat_clusters=='16')
quantile(cl16@meta.data$avg_infercnv)
cl17<-subset(nsclc.combined,subset=seurat_clusters=='17')
quantile(cl17@meta.data$avg_infercnv)
cl18<-subset(nsclc.combined,subset=seurat_clusters=='18')
quantile(cl18@meta.data$avg_infercnv)
cl19<-subset(nsclc.combined,subset=seurat_clusters=='19')
quantile(cl19@meta.data$avg_infercnv)
cl20<-subset(nsclc.combined,subset=seurat_clusters=='20')
quantile(cl20@meta.data$avg_infercnv)
cl21<-subset(nsclc.combined,subset=seurat_clusters=='21')
quantile(cl21@meta.data$avg_infercnv)
cl22<-subset(nsclc.combined,subset=seurat_clusters=='22')
quantile(cl22@meta.data$avg_infercnv)
cl23<-subset(nsclc.combined,subset=seurat_clusters=='23')
quantile(cl23@meta.data$avg_infercnv)
cl24<-subset(nsclc.combined,subset=seurat_clusters=='24')
quantile(cl24@meta.data$avg_infercnv)
cl25<-subset(nsclc.combined,subset=seurat_clusters=='25')
quantile(cl25@meta.data$avg_infercnv)
cl26<-subset(nsclc.combined,subset=seurat_clusters=='26')
quantile(cl26@meta.data$avg_infercnv)
cl27<-subset(nsclc.combined,subset=seurat_clusters=='27')
quantile(cl27@meta.data$avg_infercnv)
cl28<-subset(nsclc.combined,subset=seurat_clusters=='28')
quantile(cl28@meta.data$avg_infercnv)
cl29<-subset(nsclc.combined,subset=seurat_clusters=='29')
quantile(cl29@meta.data$avg_infercnv)
cl30<-subset(nsclc.combined,subset=seurat_clusters=='30')
quantile(cl30@meta.data$avg_infercnv)
cl31<-subset(nsclc.combined,subset=seurat_clusters=='31')
quantile(cl31@meta.data$avg_infercnv)
cl32<-subset(nsclc.combined,subset=seurat_clusters=='32')
quantile(cl32@meta.data$avg_infercnv)
cl33<-subset(nsclc.combined,subset=seurat_clusters=='33')
quantile(cl33@meta.data$avg_infercnv)
cl34<-subset(nsclc.combined,subset=seurat_clusters=='34')
quantile(cl34@meta.data$avg_infercnv)


immune.combined<-nsclc.combined

DefaultAssay(immune.combined) <- "RNA"

pdf(file = "NSCLC_KRAS_STK_samples_tumor_updated_overlay_44_samples_merge_no_anchor_test_markers_for_rare_population_infercnv_upd.pdf",height = 20,width=20)
FeaturePlot(object = immune.combined, features = 'avg_infercnv') + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

metrics <-  c("avg_infercnv")
# Extract the UMAP coordinates for each cell and include information about the metrics to plot
qc_data <- FetchData(immune.combined, vars = c(metrics, "ident", "UMAP_1", "UMAP_2"))
# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(immune.combined, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  as.data.frame() %>% group_by(ident) %>% summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plot a UMAP plot for each metric
map(metrics, function(qc){
  ggplot(qc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=qc), size=0.7) +
    scale_color_gradientn(colours = rev(rainbow(4))) +theme_bw() +
    #scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') +theme_bw() +
    geom_text(data=umap_label, aes(label=ident, x, y)) +
    ggtitle(paste0(qc,"_",pat))+
    theme(legend.text=element_text(size=7),legend.title=element_text(size=0),legend.key.width=unit(0.05,"cm"),plot.title = element_text(size=10))
}) %>% plot_grid(plotlist = .) %>% print()

dev.off()

# Average CNV: Cluster Quantile CNV
cl0<-subset(nsclc.combined,subset=seurat_clusters=='0')
quantile(cl0@meta.data$avg_infercnv)
cl1<-subset(nsclc.combined,subset=seurat_clusters=='1')
quantile(cl1@meta.data$avg_infercnv)
cl2<-subset(nsclc.combined,subset=seurat_clusters=='2')
quantile(cl2@meta.data$avg_infercnv)
cl3<-subset(nsclc.combined,subset=seurat_clusters=='3')
quantile(cl3@meta.data$avg_infercnv)
cl4<-subset(nsclc.combined,subset=seurat_clusters=='4')
quantile(cl4@meta.data$avg_infercnv)
cl5<-subset(nsclc.combined,subset=seurat_clusters=='5')
quantile(cl5@meta.data$avg_infercnv)
cl6<-subset(nsclc.combined,subset=seurat_clusters=='6')
quantile(cl6@meta.data$avg_infercnv)
cl7<-subset(nsclc.combined,subset=seurat_clusters=='7')
quantile(cl7@meta.data$avg_infercnv)
cl8<-subset(nsclc.combined,subset=seurat_clusters=='8')
quantile(cl8@meta.data$avg_infercnv)
cl9<-subset(nsclc.combined,subset=seurat_clusters=='9')
quantile(cl9@meta.data$avg_infercnv)
cl10<-subset(nsclc.combined,subset=seurat_clusters=='10')
quantile(cl10@meta.data$avg_infercnv)
cl11<-subset(nsclc.combined,subset=seurat_clusters=='11')
quantile(cl11@meta.data$avg_infercnv)
cl12<-subset(nsclc.combined,subset=seurat_clusters=='12')
quantile(cl12@meta.data$avg_infercnv)
cl13<-subset(nsclc.combined,subset=seurat_clusters=='13')
quantile(cl13@meta.data$avg_infercnv)
cl14<-subset(nsclc.combined,subset=seurat_clusters=='14')
quantile(cl14@meta.data$avg_infercnv)
cl15<-subset(nsclc.combined,subset=seurat_clusters=='15')
quantile(cl15@meta.data$avg_infercnv)
cl16<-subset(nsclc.combined,subset=seurat_clusters=='16')
quantile(cl16@meta.data$avg_infercnv)
cl17<-subset(nsclc.combined,subset=seurat_clusters=='17')
quantile(cl17@meta.data$avg_infercnv)
cl18<-subset(nsclc.combined,subset=seurat_clusters=='18')
quantile(cl18@meta.data$avg_infercnv)
cl19<-subset(nsclc.combined,subset=seurat_clusters=='19')
quantile(cl19@meta.data$avg_infercnv)
cl20<-subset(nsclc.combined,subset=seurat_clusters=='20')
quantile(cl20@meta.data$avg_infercnv)
cl21<-subset(nsclc.combined,subset=seurat_clusters=='21')
quantile(cl21@meta.data$avg_infercnv)
cl22<-subset(nsclc.combined,subset=seurat_clusters=='22')
quantile(cl22@meta.data$avg_infercnv)
cl23<-subset(nsclc.combined,subset=seurat_clusters=='23')
quantile(cl23@meta.data$avg_infercnv)
cl24<-subset(nsclc.combined,subset=seurat_clusters=='24')
quantile(cl24@meta.data$avg_infercnv)
cl25<-subset(nsclc.combined,subset=seurat_clusters=='25')
quantile(cl25@meta.data$avg_infercnv)
cl26<-subset(nsclc.combined,subset=seurat_clusters=='26')
quantile(cl26@meta.data$avg_infercnv)
cl27<-subset(nsclc.combined,subset=seurat_clusters=='27')
quantile(cl27@meta.data$avg_infercnv)
cl28<-subset(nsclc.combined,subset=seurat_clusters=='28')
quantile(cl28@meta.data$avg_infercnv)
cl29<-subset(nsclc.combined,subset=seurat_clusters=='29')
quantile(cl29@meta.data$avg_infercnv)
cl30<-subset(nsclc.combined,subset=seurat_clusters=='30')
quantile(cl30@meta.data$avg_infercnv)
cl31<-subset(nsclc.combined,subset=seurat_clusters=='31')
quantile(cl31@meta.data$avg_infercnv)
cl32<-subset(nsclc.combined,subset=seurat_clusters=='32')
quantile(cl32@meta.data$avg_infercnv)
cl33<-subset(nsclc.combined,subset=seurat_clusters=='33')
quantile(cl33@meta.data$avg_infercnv)
cl34<-subset(nsclc.combined,subset=seurat_clusters=='34')
quantile(cl34@meta.data$avg_infercnv)

# Mosaic plots
  dat <- data.frame(
    "Total_num_tumor_cells" = c(46907, 107966),
    "Cluster_21_num_tumor_cells" = c(186, 2626),
    row.names = c("Primary", "Brain_Mets"),
    stringsAsFactors = FALSE
  )
colnames(dat) <- c("Total_num_tumor_cells", "Cluster_21_num_tumor_cells")
#colnames(dat) <- c("Non-smoker", "Smoker")

dat <- data.frame(
  "Total_num_tumor_cells" = c(58572, 22113),
  "Cluster_21_num_tumor_cells" = c(255,631),
  row.names = c("STK11-mut", "Non-STK11-mut"),
  stringsAsFactors = FALSE
)

dat <- data.frame(
  "Total_num_tumor_cells" = c(58572, 74188),
  "Cluster_21_num_tumor_cells" = c(255,1926),
  row.names = c("STK11-mut", "STK11-WT"),
  stringsAsFactors = FALSE
)

dat <- data.frame(
  "Abc" = c(18, 5),
  "Cde" = c(0,15),
  row.names = c("STK11-mut", "STK11-WT"),
  stringsAsFactors = FALSE
)

dat <- data.frame(
  "Total_num_tumor_cells" = c(22113, 74188),
  "Cluster_21_num_tumor_cells" = c(631, 1926),
  row.names = c("Non-STK11-mut", "STK11-WT"),
  stringsAsFactors = FALSE
)

colnames(dat) <- c("Total_num_tumor_cells", "Cluster_21_num_tumor_cells")
#colnames(dat) <- c("Non-smoker", "Smoker")

library(grid)
pdf("Mosaic_plot_nonstkstk_vs_stkwt_cluster_21.pdf")
mosaicplot(dat,
           # gp_varnames = gpar(fontsize = 14, fontface = 1),
           # gp_labels = gpar(fontsize = 10),
           # labeling_args=list(rot_labels=c(bottom=90,top=90),gp_labels=(gpar(fontsize=14))),
           main = "Mosaic plot: Non-STK11-mut VS STK11-WT \nCluster_21 tumor cells wrt all tumor cells",
           #main = "Mosaic plot: STK11-mut VS STK11-WT \nCluster_21 tumor cells wrt all tumor cells",
           #main = "Mosaic plot: STK11-mut VS Non-STK11-mut \nCluster_21 tumor cells wrt all tumor cells",
           #main = "Mosaic plot: Primary VS BM \nCluster_21 tumor cells wrt all tumor cells",
           color = TRUE
)
dev.off()

apply(dat1, 1, 
      function(x) {
        tbl <- matrix(dat1[1:3,1:2], byrow=T)
        fisher.test(tbl, alternative="two.sided")$p.value
      })

test <- fisher.test(dat)

# Contingency tables
x <- c()
for (row in rownames(dat)) {
  for (col in colnames(dat)) {
    x <- rbind(x, matrix(rep(c(row, col), dat[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
#colnames(df) <- c("Sport_habits", "Smoking_habits")
df

# Fisher's exact test with raw data
test <- fisher.test(table(df))

# combine plot and statistical test with ggbarstats
library(ggstatsplot)
pdf("Mosaic_plot_v2_nonstk_vs_stkwt_cluster_21.pdf")
ggbarstats(
  df, V1, V2,
  results.subtitle = FALSE,
  subtitle = paste0(
    #"Mosaic plot: Primary VS BM \nCluster_21 tumor cells wrt all tumor cells\n", "Fisher's exact test", ", p-value = ",
    # "Mosaic plot: STK11-mut VS Non-STK11-mut \nCluster_21 tumor cells wrt all tumor cells\n", "Fisher's exact test", ", p-value = ",
    "Mosaic plot: Non-STK11-mut VS STK11-WT \nCluster_21 tumor cells wrt all tumor cells\n", "Fisher's exact test", ", p-value = ",
    #”Mosaic plot: STK11-mut VS STK11-WT \nCluster_21 tumor cells wrt all tumor cells\n", "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)
dev.off()

# 13) Cluster 21: DGE
# -------------------
Cluster_21_NSCLC <- readRDS(file="Cluster_21_NSCLC.rds")
colnames(Cluster_21_NSCLC@meta.data)

table(Cluster_21_NSCLC@meta.data$orig.ident)

Idents(object = NSCLC_tumor_just_merge_no_anchor) <- "seurat_clusters"
Cluster_21_0 <- FindMarkers(object = NSCLC_tumor_just_merge_no_anchor, ident.1 = 21,ident.2 = 0)
write.csv(Cluster_21_0,file="Cluster_21_0_dge.csv")
Cluster_21_1 <- FindMarkers(object = NSCLC_tumor_just_merge_no_anchor, ident.1 = 21,ident.2 = 1)
write.csv(Cluster_21_1,file="Cluster_21_1_dge.csv")
# repeat the DGE for Cluster 21 vs other clusters

# 14) Cluster 21: DGE Pairwise comparison of samples within Cluster 21
# --------------------------------------------------------------------
Idents(object = Cluster_21_NSCLC) <- "orig.ident"
PA048_KRAS_6 <- FindMarkers(object = Cluster_21_NSCLC, ident.1 = 'PA048',ident.2 = 'KRAS_6')
write.csv(PA048_KRAS_6,file="~/Documents/Ben_Izar_Project/One_on_One/12_27_2021/Rare_population_study/PA048_KRAS_6_dge.csv")
PA048_PA042 <- FindMarkers(object = Cluster_21_NSCLC, ident.1 = 'PA048',ident.2 = 'PA042')
write.csv(PA048_PA042,file="~/Documents/Ben_Izar_Project/One_on_One/12_27_2021/Rare_population_study/PA048_PA042_dge.csv")

# 15) Cluster 21: Venn
library(ggVennDiagram)
library(venn)

# List of items

cluster_21_dge <- read.csv(file="cluster_21_dge_updated_v3.1.csv")
dim(cluster_21_dge)
cluster_21_dge <- as.data.frame(cluster_21_dge)
colnames(cluster_21_dge)

x <- list('PA048_KRAS_6' = cluster_21_dge$PA048_KRAS_6,
          'PA048_PA042' = cluster_21_dge$PA048_PA042,
          'PA048_PA001' = cluster_21_dge$PA048_PA001,
          'PA048_PA068' = cluster_21_dge$PA048_PA068,
          'PA048_PA104' = cluster_21_dge$PA048_PA104,
          'PA048_STK_1' = cluster_21_dge$PA048_STK_1,
          'PA048_STK_21' = cluster_21_dge$PA048_STK_21,
          'PA048_PA067' = cluster_21_dge$PA048_PA067,
          'PA048_KRAS_13' = cluster_21_dge$PA048_KRAS_13,
          'PA048_PA080' = cluster_21_dge$PA048_PA080,
          'PA048_KRAS_7' = cluster_21_dge$PA048_KRAS_7,
          'PA048_STK_14' = cluster_21_dge$PA048_STK_14,
          'PA048_PA019' = cluster_21_dge$PA048_PA019,
          'PA048_PA034' = cluster_21_dge$PA048_PA034,
          'PA048_PA125' = cluster_21_dge$PA048_PA125,
          'PA048_KRAS_8' = cluster_21_dge$PA048_KRAS_8,
          'PA048_PA070' = cluster_21_dge$PA048_PA070,
          'PA048_KRAS_4' = cluster_21_dge$PA048_KRAS_4,
          'PA048_PA056' = cluster_21_dge$PA048_PA056,
          'PA048_PA005' = cluster_21_dge$PA048_PA005,
          'PA048_PA054' = cluster_21_dge$PA048_PA054,
          'PA048_STK_5dot2' = cluster_21_dge$PA048_STK_5dot2,
          'PA048_KRAS_12' = cluster_21_dge$PA048_KRAS_12,
          'PA048_STK_3' = cluster_21_dge$PA048_STK_3,
          'PA048_KRAS_17' = cluster_21_dge$PA048_KRAS_17,
          'PA048_N561' = cluster_21_dge$PA048_N561,
          'PA048_PA025' = cluster_21_dge$PA048_PA025,
          'PA048_PA076' = cluster_21_dge$PA048_PA076,
          'PA048_PA004' = cluster_21_dge$PA048_PA004,
          'PA048_STK_18' = cluster_21_dge$PA048_STK_18,
          'PA048_PA141' = cluster_21_dge$PA048_PA141,
          'PA048_KRAS_10' = cluster_21_dge$PA048_KRAS_10,
          'PA048_PA043' = cluster_21_dge$PA048_PA043,
          'PA048_PA072' = cluster_21_dge$PA048_PA072,
          'PA048_KRAS_11' = cluster_21_dge$PA048_KRAS_11,
          'PA048_N586' = cluster_21_dge$PA048_N586,
          'PA048_STK_15' = cluster_21_dge$PA048_STK_15,
          'PA048_STK_20' = cluster_21_dge$PA048_STK_20,
          'PA048_PA060' = cluster_21_dge$PA048_PA060,
          'KRAS_6_PA042' = cluster_21_dge$KRAS_6_PA042,
          'KRAS_6_PA001' = cluster_21_dge$KRAS_6_PA001,
          'KRAS_6_PA068' = cluster_21_dge$KRAS_6_PA068,
          'KRAS_6_PA104' = cluster_21_dge$KRAS_6_PA104,
          'KRAS_6_STK_1' = cluster_21_dge$KRAS_6_STK_1,
          'KRAS_6_STK_21' = cluster_21_dge$KRAS_6_STK_21,
          'KRAS_6_PA067' = cluster_21_dge$KRAS_6_PA067,
          'KRAS_6_KRAS_13' = cluster_21_dge$KRAS_6_KRAS_13,
          'KRAS_6_PA080' = cluster_21_dge$KRAS_6_PA080,
          'KRAS_6_KRAS_7' = cluster_21_dge$KRAS_6_KRAS_7,
          'STK_21_STK_20' = cluster_21_dge$STK_21_STK_20,
          'STK_21_STK_5dot2' = cluster_21_dge$STK_21_STK_5dot2,
          'STK_21_STK_3' = cluster_21_dge$STK_21_STK_3,
          'STK_21_PA141' = cluster_21_dge$STK_21_PA141,
          'STK_21_PA125' = cluster_21_dge$STK_21_PA125,
          'STK_21_PA080' = cluster_21_dge$STK_21_PA080,
          'STK_21_PA076' = cluster_21_dge$STK_21_PA076,
          'STK_21_PA072' = cluster_21_dge$STK_21_PA072,
          'STK_21_PA070' = cluster_21_dge$STK_21_PA070,
          'STK_21_PA067' = cluster_21_dge$STK_21_PA067,
          'STK_21_PA060' = cluster_21_dge$STK_21_PA060,
          'STK_21_PA056' = cluster_21_dge$STK_21_PA056,
          'STK_21_PA054' = cluster_21_dge$STK_21_PA054,
          'STK_21_PA043' = cluster_21_dge$STK_21_PA043,
          'STK_21_PA034' = cluster_21_dge$STK_21_PA034,
          'STK_21_PA025' = cluster_21_dge$STK_21_PA025,
          'STK_21_PA019' = cluster_21_dge$STK_21_PA019,
          'STK_21_PA005' = cluster_21_dge$STK_21_PA005,
          'STK_21_PA004' = cluster_21_dge$STK_21_PA004,
          'STK_21_N586' = cluster_21_dge$STK_21_N586,
          'STK_21_N561' = cluster_21_dge$STK_21_N561,
          'STK_21_KRAS_17' = cluster_21_dge$STK_21_KRAS_17,
          'STK_21_KRAS_13' = cluster_21_dge$STK_21_KRAS_13,
          'STK_21_KRAS_12' = cluster_21_dge$STK_21_KRAS_12,
          'STK_21_KRAS_11' = cluster_21_dge$STK_21_KRAS_11,
          'STK_21_KRAS_10' = cluster_21_dge$STK_21_KRAS_10,
          'STK_21_KRAS_8' = cluster_21_dge$STK_21_KRAS_8,
          'STK_21_KRAS_7' = cluster_21_dge$STK_21_KRAS_7,
          'STK_21_KRAS_4' = cluster_21_dge$STK_21_KRAS_4,
          'STK_20_PA060' = cluster_21_dge$STK_20_PA060,
          'STK_18_STK_20' = cluster_21_dge$STK_18_STK_20,
          'STK_18_PA141' = cluster_21_dge$STK_18_PA141,
          'STK_18_PA072' = cluster_21_dge$STK_18_PA072,
          'STK_18_PA060' = cluster_21_dge$STK_18_PA060,
          'STK_18_PA043' = cluster_21_dge$STK_18_PA043,
          'STK_18_N586' = cluster_21_dge$STK_18_N586,
          'STK_18_KRAS_11' = cluster_21_dge$STK_18_KRAS_11,
          'STK_18_KRAS_10' = cluster_21_dge$STK_18_KRAS_10,
          'STK_15_STK_20' = cluster_21_dge$STK_15_STK_20,
          'STK_15_PA060' = cluster_21_dge$STK_15_PA060,
          'STK_14_STK_20' = cluster_21_dge$STK_14_STK_20,
          'STK_14_STK_5dot2' = cluster_21_dge$STK_14_STK_5dot2,
          'STK_14_STK_3' = cluster_21_dge$STK_14_STK_3,
          'STK_14_PA141' = cluster_21_dge$STK_14_PA141,
          'STK_14_PA125' = cluster_21_dge$STK_14_PA125,
          'STK_14_PA076' = cluster_21_dge$STK_14_PA076,
          'STK_14_PA072' = cluster_21_dge$STK_14_PA072,
          'STK_14_PA070' = cluster_21_dge$STK_14_PA070,
          'STK_14_PA060' = cluster_21_dge$STK_14_PA060,
          'STK_14_PA056' = cluster_21_dge$STK_14_PA056,
          'STK_14_PA054' = cluster_21_dge$STK_14_PA054,
          'STK_14_PA043' = cluster_21_dge$STK_14_PA043,
          'STK_14_PA034' = cluster_21_dge$STK_14_PA034,
          'STK_14_PA025' = cluster_21_dge$STK_14_PA025,
          'STK_14_PA019' = cluster_21_dge$STK_14_PA019,
          'STK_14_PA005' = cluster_21_dge$STK_14_PA005,
          'STK_14_PA004' = cluster_21_dge$STK_14_PA004,
          'STK_14_N586' = cluster_21_dge$STK_14_N586,
          'STK_14_N561' = cluster_21_dge$STK_14_N561,
          'STK_14_KRAS_17' = cluster_21_dge$STK_14_KRAS_17,
          'STK_14_KRAS_12' = cluster_21_dge$STK_14_KRAS_12,
          'STK_14_KRAS_11' = cluster_21_dge$STK_14_KRAS_11,
          'STK_14_KRAS_10' = cluster_21_dge$STK_14_KRAS_10,
          'STK_14_KRAS_8' = cluster_21_dge$STK_14_KRAS_8,
          'STK_14_KRAS_4' = cluster_21_dge$STK_14_KRAS_4,
          'STK_5dot2_STK_20' = cluster_21_dge$STK_5dot2_STK_20,
          'STK_5dot2_STK_3' = cluster_21_dge$STK_5dot2_STK_3,
          'STK_5dot2_PA141' = cluster_21_dge$STK_5dot2_PA141,
          'STK_5dot2_PA076' = cluster_21_dge$STK_5dot2_PA076,
          'STK_5dot2_PA072' = cluster_21_dge$STK_5dot2_PA072,
          'STK_5dot2_PA060' = cluster_21_dge$STK_5dot2_PA060,
          'STK_5dot2_PA043' = cluster_21_dge$STK_5dot2_PA043,
          'STK_5dot2_PA025' = cluster_21_dge$STK_5dot2_PA025,
          'STK_5dot2_PA004' = cluster_21_dge$STK_5dot2_PA004,
          'STK_5dot2_N586' = cluster_21_dge$STK_5dot2_N586,
          'STK_5dot2_N561' = cluster_21_dge$STK_5dot2_N561,
          'KRAS_6_STK_14' = cluster_21_dge$KRAS_6_STK_14,
          'KRAS_6_PA019' = cluster_21_dge$KRAS_6_PA019,
          'KRAS_6_PA019' = cluster_21_dge$KRAS_6_PA019,
          'KRAS_6_PA034' = cluster_21_dge$KRAS_6_PA034,
          'KRAS_6_PA125' = cluster_21_dge$KRAS_6_PA125,
          'KRAS_6_KRAS_8' = cluster_21_dge$KRAS_6_KRAS_8,
          'KRAS_6_PA070' = cluster_21_dge$KRAS_6_PA070,
          'KRAS_6_KRAS_4' = cluster_21_dge$KRAS_6_KRAS_4,
          'KRAS_6_PA056' = cluster_21_dge$KRAS_6_PA056,
          'KRAS_6_PA005' = cluster_21_dge$KRAS_6_PA005,
          'KRAS_6_PA054' = cluster_21_dge$KRAS_6_PA054,
          'KRAS_6_STK_5dot2' = cluster_21_dge$KRAS_6_STK_5dot2,
          'KRAS_6_KRAS_12' = cluster_21_dge$KRAS_6_KRAS_12,
          'KRAS_6_STK_3' = cluster_21_dge$KRAS_6_STK_3,
          'KRAS_6_KRAS_17' = cluster_21_dge$KRAS_6_KRAS_17,
          'KRAS_6_N561' = cluster_21_dge$KRAS_6_N561,
          'KRAS_6_PA025' = cluster_21_dge$KRAS_6_PA025,
          'KRAS_6_PA076' = cluster_21_dge$KRAS_6_PA076,
          'KRAS_6_PA004' = cluster_21_dge$KRAS_6_PA004,
          'KRAS_6_STK_18' = cluster_21_dge$KRAS_6_STK_18,
          'KRAS_6_PA141' = cluster_21_dge$KRAS_6_PA141,
          'KRAS_6_KRAS_10' = cluster_21_dge$KRAS_6_KRAS_10,
          'KRAS_6_PA043' = cluster_21_dge$KRAS_6_PA043,
          'KRAS_6_PA072' = cluster_21_dge$KRAS_6_PA072,
          'KRAS_6_KRAS_11' = cluster_21_dge$KRAS_6_KRAS_11,
          'KRAS_6_N586' = cluster_21_dge$KRAS_6_N586,
          'KRAS_6_STK_15' = cluster_21_dge$KRAS_6_STK_15,
          'KRAS_6_STK_20' = cluster_21_dge$KRAS_6_STK_20,
          'KRAS_6_PA060' = cluster_21_dge$KRAS_6_PA060,
          'PA042_PA001' = cluster_21_dge$PA042_PA001,
          'PA042_PA068' = cluster_21_dge$PA042_PA068,
          'PA042_PA104' = cluster_21_dge$PA042_PA104,
          'PA042_STK_1' = cluster_21_dge$PA042_STK_1,
          'PA042_STK_21' = cluster_21_dge$PA042_STK_21,
          'PA042_PA067' = cluster_21_dge$PA042_PA067,
          'PA042_KRAS_13' = cluster_21_dge$PA042_KRAS_13,
          'PA042_PA080' = cluster_21_dge$PA042_PA080,
          'PA042_KRAS_7' = cluster_21_dge$PA042_KRAS_7,
          'PA042_STK_14' = cluster_21_dge$PA042_STK_14,
          'PA042_PA019' = cluster_21_dge$PA042_PA019,
          'PA042_PA034' = cluster_21_dge$PA042_PA034,
          'PA042_PA125' = cluster_21_dge$PA042_PA125,
          'PA042_KRAS_8' = cluster_21_dge$PA042_KRAS_8,
          'PA042_PA070' = cluster_21_dge$PA042_PA070,
          'PA042_KRAS_4' = cluster_21_dge$PA042_KRAS_4,
          'PA042_PA056' = cluster_21_dge$PA042_PA056,
          'PA042_PA005' = cluster_21_dge$PA042_PA005,
          'PA042_PA054' = cluster_21_dge$PA042_PA054,
          'PA042_STK_5dot2' = cluster_21_dge$PA042_STK_5dot2,
          'PA042_KRAS_12' = cluster_21_dge$PA042_KRAS_12,
          'PA042_STK_3' = cluster_21_dge$PA042_STK_3,
          'PA042_KRAS_17' = cluster_21_dge$PA042_KRAS_17,
          'PA042_N561' = cluster_21_dge$PA042_N561,
          'PA042_PA025' = cluster_21_dge$PA042_PA025,
          'PA042_PA076' = cluster_21_dge$PA042_PA076,
          'PA042_PA004' = cluster_21_dge$PA042_PA004,
          'PA042_STK_18' = cluster_21_dge$PA042_STK_18,
          'PA042_PA141' = cluster_21_dge$PA042_PA141,
          'PA042_KRAS_10' = cluster_21_dge$PA042_KRAS_10,
          'PA042_PA043' = cluster_21_dge$PA042_PA043,
          'PA042_PA072' = cluster_21_dge$PA042_PA072,
          'PA042_KRAS_11' = cluster_21_dge$PA042_KRAS_11,
          'PA042_N586' = cluster_21_dge$PA042_N586,
          'PA042_STK_15' = cluster_21_dge$PA042_STK_15,
          'PA042_STK_20' = cluster_21_dge$PA042_STK_20,
          'PA042_PA060' = cluster_21_dge$PA042_PA060,
          'PA001_PA068' = cluster_21_dge$PA001_PA068,
          'PA001_PA104' = cluster_21_dge$PA001_PA104,
          'PA001_STK_1' = cluster_21_dge$PA001_STK_1,
          'PA001_STK_21' = cluster_21_dge$PA001_STK_21,
          'PA001_PA067' = cluster_21_dge$PA001_PA067,
          'PA001_KRAS_13' = cluster_21_dge$PA001_KRAS_13,
          'PA001_PA080' = cluster_21_dge$PA001_PA080,
          'PA001_KRAS_7' = cluster_21_dge$PA001_KRAS_7,
          'PA001_STK_14' = cluster_21_dge$PA001_STK_14,
          'PA001_PA019' = cluster_21_dge$PA001_PA019,
          'PA001_PA034' = cluster_21_dge$PA001_PA034,
          'PA001_PA125' = cluster_21_dge$PA001_PA125,
          'PA001_KRAS_8' = cluster_21_dge$PA001_KRAS_8,
          'PA001_PA070' = cluster_21_dge$PA001_PA070,
          'PA001_KRAS_4' = cluster_21_dge$PA001_KRAS_4,
          'PA001_PA056' = cluster_21_dge$PA001_PA056,
          'PA001_PA005' = cluster_21_dge$PA001_PA005,
          'PA001_PA054' = cluster_21_dge$PA001_PA054,
          'PA001_STK_5dot2' = cluster_21_dge$PA001_STK_5dot2,
          'PA001_KRAS_12' = cluster_21_dge$PA001_KRAS_12,
          'PA001_STK_3' = cluster_21_dge$PA001_STK_3,
          'PA001_KRAS_17' = cluster_21_dge$PA001_KRAS_17,
          'PA001_N561' = cluster_21_dge$PA001_N561,
          'PA001_PA025' = cluster_21_dge$PA001_PA025,
          'PA001_PA076' = cluster_21_dge$PA001_PA076,
          'PA001_PA004' = cluster_21_dge$PA001_PA004,
          'PA001_STK_18' = cluster_21_dge$PA001_STK_18,
          'PA001_PA141' = cluster_21_dge$PA001_PA141,
          'PA001_KRAS_10' = cluster_21_dge$PA001_KRAS_10,
          'PA001_PA043' = cluster_21_dge$PA001_PA043,
          'PA001_PA072' = cluster_21_dge$PA001_PA072,
          'PA001_KRAS_11' = cluster_21_dge$PA001_KRAS_11,
          'PA001_N586' = cluster_21_dge$PA001_N586,
          'PA001_STK_15' = cluster_21_dge$PA001_STK_15,
          'PA001_STK_20' = cluster_21_dge$PA001_STK_20,
          'PA001_PA060' = cluster_21_dge$PA001_PA060,
          'PA068_PA104' = cluster_21_dge$PA068_PA104,
          'PA068_STK_1' = cluster_21_dge$PA068_STK_1,
          'PA068_STK_21' = cluster_21_dge$PA068_STK_21,
          'PA068_PA067' = cluster_21_dge$PA068_PA067,
          'PA068_KRAS_13' = cluster_21_dge$PA068_KRAS_13,
          'PA068_PA080' = cluster_21_dge$PA068_PA080,
          'PA068_KRAS_7' = cluster_21_dge$PA068_KRAS_7,
          'PA068_STK_14' = cluster_21_dge$PA068_STK_14,
          'PA068_PA019' = cluster_21_dge$PA068_PA019,
          'PA068_PA034' = cluster_21_dge$PA068_PA034,
          'PA068_PA125' = cluster_21_dge$PA068_PA125,
          'PA068_KRAS_8' = cluster_21_dge$PA068_KRAS_8,
          'PA068_PA070' = cluster_21_dge$PA068_PA070,
          'PA068_KRAS_4' = cluster_21_dge$PA068_KRAS_4,
          'PA068_PA056' = cluster_21_dge$PA068_PA056,
          'PA068_PA005' = cluster_21_dge$PA068_PA005,
          'PA068_PA054' = cluster_21_dge$PA068_PA054,
          'PA068_STK_5dot2' = cluster_21_dge$PA068_STK_5dot2,
          'PA068_KRAS_12' = cluster_21_dge$PA068_KRAS_12,
          'PA068_STK_3' = cluster_21_dge$PA068_STK_3,
          'PA068_KRAS_17' = cluster_21_dge$PA068_KRAS_17,
          'PA068_N561' = cluster_21_dge$PA068_N561,
          'PA068_PA025' = cluster_21_dge$PA068_PA025,
          'PA068_PA076' = cluster_21_dge$PA068_PA076,
          'PA068_PA004' = cluster_21_dge$PA068_PA004,
          'PA068_STK_18' = cluster_21_dge$PA068_STK_18,
          'PA068_PA141' = cluster_21_dge$PA068_PA141,
          'PA068_KRAS_10' = cluster_21_dge$PA068_KRAS_10,
          'PA068_PA043' = cluster_21_dge$PA068_PA043,
          'PA068_PA072' = cluster_21_dge$PA068_PA072,
          'PA068_KRAS_11' = cluster_21_dge$PA068_KRAS_11,
          'PA068_N586' = cluster_21_dge$PA068_N586,
          'PA068_STK_15' = cluster_21_dge$PA068_STK_15,
          'PA068_STK_20' = cluster_21_dge$PA068_STK_20,
          'PA068_PA060' = cluster_21_dge$PA068_PA060,
          'PA104_STK_1' = cluster_21_dge$PA104_STK_1,
          'PA104_STK_21' = cluster_21_dge$PA104_STK_21,
          'PA104_PA067' = cluster_21_dge$PA104_PA067,
          'PA104_KRAS_13' = cluster_21_dge$PA104_KRAS_13,
          'PA104_PA080' = cluster_21_dge$PA104_PA080,
          'PA104_KRAS_7' = cluster_21_dge$PA104_KRAS_7,
          'PA104_STK_14' = cluster_21_dge$PA104_STK_14,
          'PA104_PA019' = cluster_21_dge$PA104_PA019,
          'PA104_PA034' = cluster_21_dge$PA104_PA034,
          'PA104_PA125' = cluster_21_dge$PA104_PA125,
          'PA104_KRAS_8' = cluster_21_dge$PA104_KRAS_8,
          'PA104_PA070' = cluster_21_dge$PA104_PA070,
          'PA104_KRAS_4' = cluster_21_dge$PA104_KRAS_4,
          'PA104_PA056' = cluster_21_dge$PA104_PA056,
          'PA104_PA005' = cluster_21_dge$PA104_PA005,
          'PA104_PA054' = cluster_21_dge$PA104_PA054,
          'PA104_STK_5dot2' = cluster_21_dge$PA104_STK_5dot2,
          'PA104_KRAS_12' = cluster_21_dge$PA104_KRAS_12,
          'PA104_STK_3' = cluster_21_dge$PA104_STK_3,
          'PA104_KRAS_17' = cluster_21_dge$PA104_KRAS_17,
          'PA104_N561' = cluster_21_dge$PA104_N561,
          'PA104_PA025' = cluster_21_dge$PA104_PA025,
          'PA104_PA076' = cluster_21_dge$PA104_PA076,
          'PA104_PA004' = cluster_21_dge$PA104_PA004,
          'PA104_STK_18' = cluster_21_dge$PA104_STK_18,
          'PA104_PA141' = cluster_21_dge$PA104_PA141,
          'PA104_KRAS_10' = cluster_21_dge$PA104_KRAS_10,
          'PA104_PA043' = cluster_21_dge$PA104_PA043,
          'PA104_PA072' = cluster_21_dge$PA104_PA072,
          'PA104_KRAS_11' = cluster_21_dge$PA104_KRAS_11,
          'PA104_N586' = cluster_21_dge$PA104_N586,
          'PA104_STK_15' = cluster_21_dge$PA104_STK_15,
          'PA104_STK_20' = cluster_21_dge$PA104_STK_20,
          'PA104_PA060' = cluster_21_dge$PA104_PA060,
          'STK_1_STK_21' = cluster_21_dge$STK_1_STK_21,
          'STK_1_PA067' = cluster_21_dge$STK_1_PA067,
          'STK_1_KRAS_13' = cluster_21_dge$STK_1_KRAS_13,
          'STK_1_PA080' = cluster_21_dge$STK_1_PA080,
          'STK_1_KRAS_7' = cluster_21_dge$STK_1_KRAS_7,
          'STK_1_STK_14' = cluster_21_dge$STK_1_STK_14,
          'STK_1_PA019' = cluster_21_dge$STK_1_PA019,
          'STK_1_PA034' = cluster_21_dge$STK_1_PA034,
          'STK_1_PA125' = cluster_21_dge$STK_1_PA125,
          'STK_1_KRAS_8' = cluster_21_dge$STK_1_KRAS_8,
          'STK_1_PA070' = cluster_21_dge$STK_1_PA070,
          'STK_1_KRAS_4' = cluster_21_dge$STK_1_KRAS_4,
          'STK_1_PA056' = cluster_21_dge$STK_1_PA056,
          'STK_1_PA005' = cluster_21_dge$STK_1_PA005,
          'STK_1_PA054' = cluster_21_dge$STK_1_PA054,
          'STK_1_STK_5dot2' = cluster_21_dge$STK_1_STK_5dot2,
          'STK_1_KRAS_12' = cluster_21_dge$STK_1_KRAS_12,
          'STK_1_STK_3' = cluster_21_dge$STK_1_STK_3,
          'STK_1_KRAS_17' = cluster_21_dge$STK_1_KRAS_17,
          'STK_1_N561' = cluster_21_dge$STK_1_N561,
          'STK_1_PA025' = cluster_21_dge$STK_1_PA025,
          'STK_1_PA076' = cluster_21_dge$STK_1_PA076,
          'STK_1_PA004' = cluster_21_dge$STK_1_PA004,
          'STK_1_STK_18' = cluster_21_dge$STK_1_STK_18,
          'STK_1_PA141' = cluster_21_dge$STK_1_PA141,
          'STK_1_KRAS_10' = cluster_21_dge$STK_1_KRAS_10,
          'STK_1_PA043' = cluster_21_dge$STK_1_PA043,
          'STK_1_PA072' = cluster_21_dge$STK_1_PA072,
          'STK_1_KRAS_11' = cluster_21_dge$STK_1_KRAS_11,
          'STK_1_N586' = cluster_21_dge$STK_1_N586,
          'STK_1_STK_15' = cluster_21_dge$STK_1_STK_15,
          'STK_1_STK_20' = cluster_21_dge$STK_1_STK_20,
          'STK_1_PA060' = cluster_21_dge$STK_1_PA060,
          'PA067_KRAS_13' = cluster_21_dge$PA067_KRAS_13,
          'PA067_PA080' = cluster_21_dge$PA067_PA080,
          'PA067_KRAS_7' = cluster_21_dge$PA067_KRAS_7,
          'PA067_PA019' = cluster_21_dge$PA067_PA019,
          'PA067_PA034' = cluster_21_dge$PA067_PA034,
          'PA067_PA125' = cluster_21_dge$PA067_PA125,
          'PA067_KRAS_8' = cluster_21_dge$PA067_KRAS_8,
          'PA067_PA070' = cluster_21_dge$PA067_PA070,
          'PA067_KRAS_4' = cluster_21_dge$PA067_KRAS_4,
          'PA067_PA056' = cluster_21_dge$PA067_PA056,
          'PA067_PA005' = cluster_21_dge$PA067_PA005,
          'PA067_PA054' = cluster_21_dge$PA067_PA054,
          'PA067_STK_5dot2' = cluster_21_dge$PA067_STK_5dot2,
          'PA067_KRAS_12' = cluster_21_dge$PA067_KRAS_12,
          'PA067_STK_3' = cluster_21_dge$PA067_STK_3,
          'PA067_KRAS_17' = cluster_21_dge$PA067_KRAS_17,
          'PA067_N561' = cluster_21_dge$PA067_N561,
          'PA067_PA025' = cluster_21_dge$PA067_PA025,
          'PA067_PA076' = cluster_21_dge$PA067_PA076,
          'PA067_PA004' = cluster_21_dge$PA067_PA004,
          'PA067_PA141' = cluster_21_dge$PA067_PA141,
          'PA067_KRAS_10' = cluster_21_dge$PA067_KRAS_10,
          'PA067_PA043' = cluster_21_dge$PA067_PA043,
          'PA067_PA072' = cluster_21_dge$PA067_PA072,
          'PA067_KRAS_11' = cluster_21_dge$PA067_KRAS_11,
          'PA067_N586' = cluster_21_dge$PA067_N586,
          'PA067_STK_20' = cluster_21_dge$PA067_STK_20,
          'PA067_PA060' = cluster_21_dge$PA067_PA060,
          'KRAS_13_PA080' = cluster_21_dge$KRAS_13_PA080,
          'KRAS_13_KRAS_7' = cluster_21_dge$KRAS_13_KRAS_7,
          'KRAS_13_PA019' = cluster_21_dge$KRAS_13_PA019,
          'KRAS_13_PA034' = cluster_21_dge$KRAS_13_PA034,
          'KRAS_13_PA125' = cluster_21_dge$KRAS_13_PA125,
          'KRAS_13_KRAS_8' = cluster_21_dge$KRAS_13_KRAS_8,
          'KRAS_13_PA070' = cluster_21_dge$KRAS_13_PA070,
          'KRAS_13_KRAS_4' = cluster_21_dge$KRAS_13_KRAS_4,
          'KRAS_13_PA056' = cluster_21_dge$KRAS_13_PA056,
          'N586_PA060' = cluster_21_dge$N586_PA060,
          'N586_STK_20' = cluster_21_dge$N586_STK_20,
          'N586_STK_15' = cluster_21_dge$N586_STK_15,
          'KRAS_11_PA060' = cluster_21_dge$KRAS_11_PA060,
          'KRAS_11_STK_20' = cluster_21_dge$KRAS_11_STK_20,
          'KRAS_11_N586' = cluster_21_dge$KRAS_11_N586,
          'PA072_PA060' = cluster_21_dge$PA072_PA060,
          'PA072_STK_20' = cluster_21_dge$PA072_STK_20,
          'PA072_N586' = cluster_21_dge$PA072_N586,
          'PA072_KRAS_11' = cluster_21_dge$PA072_KRAS_11,
          'PA043_PA060' = cluster_21_dge$PA043_PA060,
          'PA043_STK_20' = cluster_21_dge$PA043_STK_20,
          'PA043_N586' = cluster_21_dge$PA043_N586,
          'PA043_KRAS_11' = cluster_21_dge$PA043_KRAS_11,
          'PA043_PA072' = cluster_21_dge$PA043_PA072,
          'KRAS_10_PA060' = cluster_21_dge$KRAS_10_PA060,
          'KRAS_10_STK_20' = cluster_21_dge$KRAS_10_STK_20,
          'KRAS_10_N586' = cluster_21_dge$KRAS_10_N586,
          'KRAS_10_KRAS_11' = cluster_21_dge$KRAS_10_KRAS_11,
          'KRAS_10_PA072' = cluster_21_dge$KRAS_10_PA072,
          'KRAS_10_PA043' = cluster_21_dge$KRAS_10_PA043,
          'PA141_PA060' = cluster_21_dge$PA141_PA060,
          'PA141_STK_20' = cluster_21_dge$PA141_STK_20,
          'PA141_N586' = cluster_21_dge$PA141_N586,
          'PA141_KRAS_11' = cluster_21_dge$PA141_KRAS_11,
          'PA141_PA072' = cluster_21_dge$PA141_PA072,
          'PA141_PA043' = cluster_21_dge$PA141_PA043,
          'PA141_KRAS_10' = cluster_21_dge$PA141_KRAS_10,
          'PA004_PA060' = cluster_21_dge$PA004_PA060,
          'PA004_STK_20' = cluster_21_dge$PA004_STK_20,
          'PA004_N586' = cluster_21_dge$PA004_N586,
          'PA004_KRAS_11' = cluster_21_dge$PA004_KRAS_11,
          'PA004_PA072' = cluster_21_dge$PA004_PA072,
          'PA004_PA043' = cluster_21_dge$PA004_PA043,
          'PA004_KRAS_10' = cluster_21_dge$PA004_KRAS_10,
          'PA004_PA141' = cluster_21_dge$PA004_PA141,
          'PA004_STK_18' = cluster_21_dge$PA004_STK_18,
          'PA076_PA060' = cluster_21_dge$PA076_PA060,
          'PA076_STK_20' = cluster_21_dge$PA076_STK_20,
          'PA076_N586' = cluster_21_dge$PA076_N586,
          'PA076_KRAS_11' = cluster_21_dge$PA076_KRAS_11,
          'PA076_PA072' = cluster_21_dge$PA076_PA072,
          'PA076_PA043' = cluster_21_dge$PA076_PA043,
          'PA076_KRAS_10' = cluster_21_dge$PA076_KRAS_10,
          'PA076_PA141' = cluster_21_dge$PA076_PA141,
          'PA076_STK_18' = cluster_21_dge$PA076_STK_18,
          'PA076_PA004' = cluster_21_dge$PA076_PA004,
          'PA025_PA060' = cluster_21_dge$PA025_PA060,
          'PA025_STK_20' = cluster_21_dge$PA025_STK_20,
          'PA025_N586' = cluster_21_dge$PA025_N586,
          'PA025_KRAS_11' = cluster_21_dge$PA025_KRAS_11,
          'PA025_PA072' = cluster_21_dge$PA025_PA072,
          'PA025_PA043' = cluster_21_dge$PA025_PA043,
          'PA025_KRAS_10' = cluster_21_dge$PA025_KRAS_10,
          'PA025_PA141' = cluster_21_dge$PA025_PA141,
          'PA025_PA004' = cluster_21_dge$PA025_PA004,
          'PA025_PA076' = cluster_21_dge$PA025_PA076,
          'N561_PA060' = cluster_21_dge$N561_PA060,
          'N561_STK_20' = cluster_21_dge$N561_STK_20,
          'N561_N586' = cluster_21_dge$N561_N586,
          'N561_KRAS_11' = cluster_21_dge$N561_KRAS_11,
          'N561_PA072' = cluster_21_dge$N561_PA072,
          'N561_PA043' = cluster_21_dge$N561_PA043,
          'N561_KRAS_10' = cluster_21_dge$N561_KRAS_10,
          'N561_PA141' = cluster_21_dge$N561_PA141,
          'N561_PA004' = cluster_21_dge$N561_PA004,
          'N561_PA076' = cluster_21_dge$N561_PA076,
          'N561_PA025' = cluster_21_dge$N561_PA025,
          'KRAS_13_PA005' = cluster_21_dge$KRAS_13_PA005,
          'KRAS_13_PA054' = cluster_21_dge$KRAS_13_PA054,
          'PA080_KRAS_7' = cluster_21_dge$PA080_KRAS_7,
          'PA080_STK_14' = cluster_21_dge$PA080_STK_14,
          'KRAS_7_STK_14' = cluster_21_dge$KRAS_7_STK_14,
          'KRAS_7_PA070' = cluster_21_dge$KRAS_7_PA070,
          'PA019_PA034' = cluster_21_dge$PA019_PA034,
          'PA019_PA056' = cluster_21_dge$PA019_PA056,
          'PA034_PA125' = cluster_21_dge$PA034_PA125,
          'PA034_KRAS_12' = cluster_21_dge$PA034_KRAS_12,
          'PA125_PA005' = cluster_21_dge$PA125_PA005,
          'PA125_PA076' = cluster_21_dge$PA125_PA076,
          'KRAS_8_PA056' = cluster_21_dge$KRAS_8_PA056,
          'KRAS_8_N561' = cluster_21_dge$KRAS_8_N561,
          'PA070_N561' = cluster_21_dge$PA070_N561,
          'PA070_KRAS_11' = cluster_21_dge$PA070_KRAS_11,
          'KRAS_4_STK_3' = cluster_21_dge$KRAS_4_STK_3,
          'KRAS_4_PA004' = cluster_21_dge$KRAS_4_PA004,
          'PA056_PA054' = cluster_21_dge$PA056_PA054,
          'PA056_PA076' = cluster_21_dge$PA056_PA076,
          'PA005_PA025' = cluster_21_dge$PA005_PA025,
          'PA005_PA072' = cluster_21_dge$PA005_PA072,
          'PA054_PA025' = cluster_21_dge$PA054_PA025,
          'PA054_PA043' = cluster_21_dge$PA054_PA043,
          'STK_5dot2_KRAS_12' = cluster_21_dge$STK_5dot2_KRAS_12,
          'STK_5dot2_KRAS_11' = cluster_21_dge$STK_5dot2_KRAS_11,
          'KRAS_12_PA025' = cluster_21_dge$KRAS_12_PA025,
          'KRAS_12_PA060' = cluster_21_dge$KRAS_12_PA060,
          'STK_3_PA025' = cluster_21_dge$STK_3_PA025,
          'STK_3_KRAS_11' = cluster_21_dge$STK_3_KRAS_11,
          'KRAS_17_PA025' = cluster_21_dge$KRAS_17_PA025,
          'KRAS_17_PA141' = cluster_21_dge$KRAS_17_PA141,
          'KRAS_13_PA025' = cluster_21_dge$KRAS_13_PA025,
          'PA080_PA019' = cluster_21_dge$PA080_PA019,
          'KRAS_7_PA019' = cluster_21_dge$KRAS_7_PA019,
          'PA019_PA125' = cluster_21_dge$PA019_PA125,
          'PA034_PA070' = cluster_21_dge$PA034_PA070,
          'PA125_PA070' = cluster_21_dge$PA125_PA070,
          'KRAS_8_PA054' = cluster_21_dge$KRAS_8_PA054,
          'PA070_KRAS_4' = cluster_21_dge$PA070_KRAS_4,
          'KRAS_4_PA005' = cluster_21_dge$KRAS_4_PA005,
          'PA056_PA005' = cluster_21_dge$PA056_PA005,
          'PA005_PA054' = cluster_21_dge$PA005_PA054,
          'PA054_PA025' = cluster_21_dge$PA054_PA025,
          'STK_5dot2_KRAS_10' = cluster_21_dge$STK_5dot2_KRAS_10,
          'KRAS_12_PA076' = cluster_21_dge$KRAS_12_PA076,
          'STK_3_PA076' = cluster_21_dge$STK_3_PA076,
          'KRAS_17_N561' = cluster_21_dge$KRAS_17_N561,
          'KRAS_13_KRAS_12' = cluster_21_dge$KRAS_13_KRAS_12,
          'PA080_PA034' = cluster_21_dge$PA080_PA034,
          'KRAS_7_PA034' = cluster_21_dge$KRAS_7_PA034,
          'PA019_PA070' = cluster_21_dge$PA019_PA070,
          'PA034_KRAS_4' = cluster_21_dge$PA034_KRAS_4,
          'PA125_PA056' = cluster_21_dge$PA125_PA056,
          'PA125_KRAS_12' = cluster_21_dge$PA125_KRAS_12,
          'KRAS_8_PA005' = cluster_21_dge$KRAS_8_PA005,
          'KRAS_8_KRAS_11' = cluster_21_dge$KRAS_8_KRAS_11,
          'PA070_PA054' = cluster_21_dge$PA070_PA054,
          'KRAS_4_PA056' = cluster_21_dge$KRAS_4_PA056,
          'PA056_KRAS_12' = cluster_21_dge$PA056_KRAS_12,
          'PA005_STK_3' = cluster_21_dge$PA005_STK_3,
          'PA054_PA076' = cluster_21_dge$PA054_PA076,
          'STK_5dot2_KRAS_17' = cluster_21_dge$STK_5dot2_KRAS_17,
          'KRAS_12_PA004' = cluster_21_dge$KRAS_12_PA004,
          'STK_3_PA004' = cluster_21_dge$STK_3_PA004,
          'KRAS_17_PA043' = cluster_21_dge$KRAS_17_PA043,
          'PA080_PA125' = cluster_21_dge$PA080_PA125,
          'PA019_PA054' = cluster_21_dge$PA019_PA054)
# dim(x)
x
# Venn diagram with count labels
ggVennDiagram(x, 
              label = "count")

venn(x)
library(UpSetR)
# #install.packages("UpSetR")
# 
# 
# upset(fromExpression(x), 
#       nintersects = 40, 
#       nsets = 6, 
#       order.by = "freq", 
#       decreasing = T, 
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0, 
#       text.scale = 1.1, 
#       point.size = 2.8, 
#       line.size = 1
# )
write.csv(fromList(x),file="fromList.x.csv")
pdf("cluster_21_upset_top150.pdf",height = 45,width = 20)
upset(fromList(x),nsets = 150, nintersects = 100,order.by = "freq")
dev.off()

library("SuperExactTest")
res=supertest(x, n=20000)

plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15))

# 16) Cluster 21: Volcano plot

library(ggplot2)
library(scales)

cluster_21<-read.csv(file="cluster_21_dge.csv")
dim(cluster_21)
cluster_21[1:15,]
dim(tmp)
tmp[1:5,]

# remove rows that contain NA values
de <- cluster_21[complete.cases(cluster_21), ]
dim(de)


# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"

library(EnhancedVolcano)

# define different cell-types that will be shaded
celltype1 <- c('AREG','EPCAM','GFAP','HOPX','LRRK2','MBP','NKX2-1','NRG3','PCDH9','PLP1','PLXDC2','SFTA3',
               'SFTPB','SPARCL1','TCF4','VIM','ZEB1','ZEB2')
#celltype2 <- c('SORT1', 'KLF15')
has_ggalt <- ! is(try(find.package("ggalt")), "try-error")

keyvals <- ifelse(
  de$log2FoldChange < -1.0, 'royalblue',
  ifelse(de$log2FoldChange > 1.0, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'
pdf(file="volcano_cluster_21_upd.pdf",height=10,width = 10)
EnhancedVolcano(de,
                lab = de$gene_symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(if (has_ggalt) celltype1 else NULL),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Cluster_21: Rare population',
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                shape = 16,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 20,
                legendIconSize = 20.0,
                # encircle
                encircle = if (has_ggalt) celltype1 else NULL,
                encircleCol = 'black',
                encircleSize = 2.5,
                encircleFill = 'pink',
                encircleAlpha = 1/2,
                #pointSize=2,
                # shade
                #shade = celltype2,
                shadeAlpha = 1/2,
                shadeFill = 'skyblue',
                shadeSize = 1,
                shadeBins = 5,
                #  transcriptLabSize = 5,
                #drawConnectors = TRUE,
                widthConnectors = 0.2,
                drawConnectors = TRUE,
                #widthConnectors = 2.0,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 5,
                borderColour = 'black')
dev.off()

