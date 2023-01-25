## ---------------------------
##
## Purpose of script: Calculates DEGs for Brain Metastasis vs Lung, STK mutant vs NonSTK Mutant
##
## Author: Abhi Jaiswal
##
## Date Created: 2022-05-30
##
## Email: ajaiswal1995@gmail.com
##
## Packages:

library(Seurat)
library(msigdbi)

####################################################################################################################################
## THIS SECTION WAS RAN ON CANNES - SEURAT VERSION 4.1.1

TCellObject = readRDS("../NSCLC_44_samples_T_cells_v3.rds")


TCellObject_CD8 = subset(TCellObject,cells = WhichCells(TCellObject, ident = c("CD8+"), slot = "manual.annot"))
TCellObject_CD8 = subset(TCellObject_CD8,cells = rownames(TCellObject_CD8@meta.data)[!TCellObject_CD8@meta.data$orig.ident %in% "N561"])


DEGs_CD8_BMvsPrimary = FindMarkers(TCellObject_CD8, ident.1 = "BRAIN_METS",ident.2 = "PRIMARY",
                                      group.by = "prim_vs_bm", min.pct = 0.25)

TCellObject_LocationSplit = SplitObject(TCellObject_CD8, split.by = "prim_vs_bm")

DEGs_CD8_STKvsNonSTKMut = lapply(TCellObject_LocationSplit,function(x) {
  FindMarkers(x, ident.1 = "STK11-mut",ident.2 = "Non-STK11-mut",
              group.by = "stk_vs_nonstk", min.pct = 0.25)
})

save(DEGs_CD8_BMvsPrimary,DEGs_CD8_STKvsNonSTKMut, file = "Izar_CD8TCells_DEGAnalysis_v3.rda")

####################################################################################################################################
## THIS SECTION WAS RAN LOCALLY

load("Izar_CD8TCells_DEGAnalysis_v3.rda")

DEGList = list(DEGs_CD8_BMvsPrimary = DEGs_CD8_BMvsPrimary,
               DEGs_CD8_Primary_STKvsNonSTKMut = DEGs_CD8_STKvsNonSTKMut$PRIMARY,
               DEGs_CD8_BM_STKvsNonSTKMut = DEGs_CD8_STKvsNonSTKMut$BRAIN_METS)

FDRCutoff = 0.05
LogFCCutoff = 0

DEGList = lapply(DEGList, function(x) {
  
  Up = rownames(subset(x, avg_log2FC > LogFCCutoff & p_val_adj < FDRCutoff))
  
  Dn = rownames(subset(x, avg_log2FC < (-1*LogFCCutoff) & p_val_adj < FDRCutoff))
  
  list(Up = Up, Dn = Dn)
})

DEGList = unlist(DEGList,recursive = F)

DEGList = list(genesets = DEGList,
               geneset.names = names(DEGList),
               geneset.descriptions = names(DEGList))

msigdbi::write.gmt(DEGList,"Izar_CD8TCell_DEGs_v3.gmt")
