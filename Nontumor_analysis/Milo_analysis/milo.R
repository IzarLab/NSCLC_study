# Author: Somnath Tagore, Ph.D. Title: Milo analysis
# Script Name: milo.R 
# Last Updated: 06/24/2022

# Packages required for this analysis
#BiocManager::install("miloR")
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)
library(beeswarm)

library(rstatix)
library(ggplot2)
#install.packages("ggprism")
library(ggprism)
library(patchwork)
library(magrittr)
library(ggpubr)

# Read data
pbmc_small <- readRDS(file="./NSCLC_44_samples_T_cells_v6.rds")
colnames(pbmc_small@meta.data)
pbmc_small[['barcodes']] <- rownames(pbmc_small@meta.data)

# Create a single cell experiment object
pbmc_small_sce <- as.SingleCellExperiment(pbmc_small)
pbmc_small_sce

# Run Milo for test
pbmc_small_milo <- Milo(pbmc_small_sce)

# build graph
traj_milo <- buildGraph(pbmc_small_milo, k = 50, d = 30)

# create neighborhoods
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 50, d=30, refined = TRUE)

# plot neigborhood size
pdf(file="t_cells_milo_plotNhoodSizeHist_k50d30.pdf")
plotNhoodSizeHist(traj_milo)
dev.off()

# trajectory design
#traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample='manual.annot.fine')
#traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample='orig.ident')
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample='orig.ident')
#head(nhoodCounts(traj_milo))
traj_milo

traj_design <- data.frame(colData(traj_milo))[,c("orig.ident", "manual.annot.fine", "PRIMARY_vs_BRAIN_METS", "STK11mut_vs_NonSTK11mut","ident")]
#traj_design$manual.annot.fine <- as.factor(traj_design$manual.annot.fine) 
traj_design$manual.annot.fine <- as.factor(traj_design$manual.annot.fine) 
traj_design <- distinct(traj_design)

write.csv(traj_design,file="./traj_design_k50d30_upd.csv")

#embryo_design and test neighborhoods
traj_design <- distinct(traj_design)

traj_milo <- calcNhoodDistance(traj_milo, d=30)
#rownames(traj_design)=traj_design[,"manual.annot.fine"]

#da_results <- testNhoods(traj_milo, design = ~ manual.annot.fine, design.df = traj_design)
da_results <- testNhoods(traj_milo, design = ~ PRIMARY_vs_BRAIN_METS, design.df = traj_design)

write.csv(traj_design,file="traj_design_tcells_k50d30.csv")
dim(traj_milo@nhoodCounts)
traj_milo@nhoodCounts[1:5,1:5]

# arrange data as per spatial FDR
da_results %>%
  arrange(- SpatialFDR) %>%
  head()

pdf(file="t_cells_milo_da_results_histogram_k50d30.pdf")
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 0.03) 
dev.off()

# annotate neighborhoods as per manual annotation
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "manual.annot.fine")
#da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "manual.annot")
#da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "cell_type_fine")
head(da_results)

pdf(file="t_cells_milo_da_results_single_cell_UMAP_k10d30.pdf")
umap_pl <- plotReducedDim(traj_milo, dimred = "UMAP", colour_by="manual.annot.fine", text_by = "manual.annot.fine", 
 #                         umap_pl <- plotReducedDim(traj_milo, dimred = "UMAP", colour_by="manual.annot", text_by = "manual.annot", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
dev.off()

# build neighborhood graph
traj_milo <- buildNhoodGraph(traj_milo)

## Plot neighbourhood graph
pdf(file="t_cells_milo_da_results_neighbourhood_graph_k50d30.pdf")
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results, layout="UMAP",alpha=0.05) 
  
umap_pl + nh_graph_pl + plot_layout(guides="collect")
dev.off()

# Optional
da_results$manual.annot.fine <- as.factor(da_results$manual.annot.fine)
#da_results$cell_type_fine <- as.factor(da_results$cell_type_fine)
da_results$manual.annot.fine_fraction_1 <- ifelse(da_results$manual.annot.fine_fraction < 0.5, "Not_significant", "Significant")
da_results$manual.annot.fine_fraction_color <- ifelse(da_results$manual.annot.fine_fraction < 0.5, "Red", "Blue")

# Bee swarm plot by group
pdf(file="t_cells_milo_da_results_plotDAbeeswarm_k50d30.pdf",height=10,width=40)
beeswarm(logFC ~ manual.annot.fine, data = da_results,
    method = 'swarm',
    pch = 16, pwcol = as.numeric(manual.annot.fine_fraction_1),
    xlab = '', ylab = 'Follow-up time (months)',
    labels = c('Censored', 'Metastasis'))
  legend('topright', legend = levels(da_results$manual.annot.fine_fraction_1),
    title = 'ER', pch = 16, col = 1:2)

beeswarm(logFC ~ manual.annot.fine, data=da_results,
         pch = 19,
         pwcol = manual.annot.fine_fraction_color,horizontal=FALSE)

dev.off()

## Plot neighbourhood graph
pdf(file="t_cells_milo_da_results_plotDAbeeswarm_k50d30.pdf")
plotDAbeeswarm(da_results, group.by = "manual.annot.fine_fraction")
dev.off()

# Finer neighborhoods
da_results$celltype <- ifelse(da_results$manual.annot.fine_fraction < 0.7, "Mixed", da_results$manual.annot.fine_fraction)

da_results <- groupNhoods(traj_milo, da_results)
head(da_results)
write.csv(da_results, file="t_Cells_milo_fractions_k_100.csv")

traj_milo <- buildNhoodGraph(traj_milo)
pdf(file="t_cells_milo_1.1.pdf")

plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

dev.off()

pdf(file="t_cells_milo_3.pdf")
#plotNhoodGroups(traj_milo, da_results, layout="umap") 
plotDAbeeswarm(da_results, group.by = "manual.annot.fine")
#plotNhoodGroups(traj_milo, da_results, layout="umap") 
dev.off()

# milo fractions

tumor_nontumor <- read.csv(file="~/Documents/Ben_Izar_Project/milo/t-cells/t_cells_prim_vs_bm.csv")

row.names(tumor_nontumor) <- tumor_nontumor[,1]
tumor_nontumor<-tumor_nontumor[,-1]
rownames(tumor_nontumor)
tumor_nontumor<-as.data.frame(tumor_nontumor)

#color codes
cols = c(
  'B_cell' ='#808080',
  'Macrophage_M1' ='#d3d3d3',
  'Macrophage_M2' ='#2f4f4f',
  'Monocyte' ='#556b2f',
  'Neutrophil' ='#8b4513',
  'NK_cell' ='#2e8b57',
  'T_cell_CD4+_(non_regulatory)' ='#191970',
  'T_cell_CD8+' ='#87cefa',
  'T_cell_regulatory_(Tregs)' ='#dc143c',
  'Myeloid_dendritic_cell' ='#6a5acd',
  'uncharacterized_cell' ='#32cd32'
)


 cols = c(
  'PRIMARY' ='red',
  'BRAIN_METS' ='blue')

cols = c(
  'STK11-MUT' ='brown',
  'STK11-WT' ='grey')

tumor_nontumor$sdv <- c(0.303144704,
                        0.309818139,
                        0.317108257,
                        0.312389479,
                        0.319343209,
                        0.327656029,
                        0.33425399,
                        0.341070669,
                        0.340372989,
                        0.337044305,
                        0.330011423,
                        0.315933947)
p<-ggplot(data = tumor_nontumor, aes(x = Celltype, y = Milo_Fraction, fill = PRIMARY_VS_BM)) + 
  geom_bar( stat='identity', position=position_dodge())+theme_bw()+scale_fill_manual(values=cols)+geom_text(
  aes(label = round(Milo_Fraction,2), y =round(Milo_Fraction,2) + 0.05),
  position = position_dodge(0.9),
  vjust = 0.25
) +ggtitle("Milo Fractions (T_cells): STK11-MUT vs STK11-WT")+ #ggtitle("Milo Fractions (T_cells): PRIMARY vs BM")
   xlab("Cell types") + ylab("Milo Fractions (percentage)")+theme(axis.text.x=element_text(size=rel(1.1)))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

 ggsave("./Milo_Fractions_T_cells_PRIMARY_BM.pdf",height = 10,width = 25)
 

 #accuracy
  cols1 = c(
   'N' ='#808080',
   'Sig_N' ='#d3d3d3',
   'NonSig_N' ='#2f4f4f',
   'Accuracy' ='#556b2f')
   
 tumor_nontumor <- read.csv(file="./t_cells_milo_accuracy.csv")
 #tumor_nontumor <- read.csv(file="STK11mut_vs_STK11wt_LUAD_TCGA_deconvolution_v2.csv")
 row.names(tumor_nontumor) <- tumor_nontumor[,1]
 tumor_nontumor<-tumor_nontumor[,-1]
 rownames(tumor_nontumor)
 tumor_nontumor<-as.data.frame(tumor_nontumor)

 ggplot(data = tumor_nontumor, aes(x = Celltype, y = Prediction, fill = Features)) + 
   geom_bar( stat='identity', position=position_dodge())+theme_bw()+scale_fill_manual(values=cols1)+geom_text(
     aes(label = round(Prediction,2), y =round(Prediction,2) + 0.9),
     position = position_dodge(0.9),
     vjust = 0.25
   ) +ggtitle("Milo Predictions: T cells")+ 
     xlab("Cell types") + ylab("Milo Predictions")+theme(axis.text.x=element_text(size=rel(1.1)))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 ggsave("~/Documents/Ben_Izar_Project/milo/t-cells/Milo_Fractions_T_cells_N.pdf",height = 10,width = 25)
 
 
 # beeswarm /stripcharts final
 
 tumor_nontumor <- read.csv(file="./t_cells_milo_da_results_upda_k50d30_fraction_bm.stkmut.vs.stkwt.csv")
 row.names(tumor_nontumor) <- tumor_nontumor[,1]
 tumor_nontumor<-tumor_nontumor[,-1]
 rownames(tumor_nontumor)
 tumor_nontumor<-as.data.frame(tumor_nontumor)
 
 colnames(tumor_nontumor)
 #N = 100
 DF = tibble(Celltype = tumor_nontumor$Celltype,
             logFC = tumor_nontumor$logFC,
             SpatialFDR = tumor_nontumor$SpatialFDR)#,
 
 gg1<-ggplot(DF, aes(logFC, Celltype, fill = Celltype)) +
   geom_jitter(aes(size = SpatialFDR), shape = 21, alpha = .8, width = .05, height = .05, stroke = 0) +
   #geom_jitter(aes(size = Neighborhood_ratio_median),shape = 21, alpha = .8, width = .05, height = .05, stroke = 0) +
   theme_minimal() +
  theme(legend.position = "right")+ggtitle("Milo Fractions (T-cells): Primary vs BM")
 
 
