library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)
library(celldex)
library(SingleR)

#options(future.globals.maxSize= 50504148480)

#filt.list <- readRDS('ObjectsToIntegrate.rds')
#features <- SelectIntegrationFeatures(object.list = filt.list, nfeatures = 3000)

#filt.list <- lapply(X = filt.list, FUN = RunPCA, features = features, npcs = 40)


#immune.anchors <- FindIntegrationAnchors(object.list = filt.list, normalization.method = "SCT",anchor.features = features, dims = 1:30, reduction = "rpca")

#saveRDS(immune.anchors, "AnchorsTumorCtrlIntegratedAssay.rds")          

immune.anchors <-readRDS('AnchorsTumorCtrlIntegratedAssay.rds')

immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30, k.weight = 40)

saveRDS(immune.combined.sct, "FullMergedTumorCtrlIntegratedAssay.rds")


immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
#saveRDS(immune.combined, "FullMergedTumorCtrlIntegratedAssay.rds")
