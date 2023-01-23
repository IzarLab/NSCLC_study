## Conversion Script

library(Seurat)
library(SeuratData)
library(SeuratDisk)

rds_temp <- readRDS("data/raw_data/Control_lung_samples/data_C51ctr_cb.rds")
SaveH5Seurat(rds_temp, filename = "data/Seurat2AdataIntermediateFiles/C51_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C51_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = rds_temp, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c51.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

c52 <- readRDS("data/raw_data/Control_lung_samples/data_C52ctr_cb.rds")
SaveH5Seurat(c52, filename = "data/Seurat2AdataIntermediateFiles/C52_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C52_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c52, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c52.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


c53 <- readRDS("data/raw_data/Control_lung_samples/data_C53ctr_cb.rds")
SaveH5Seurat(c53, filename = "data/Seurat2AdataIntermediateFiles/C53_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C53_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c53, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c53.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


c54 <- readRDS("data/raw_data/Control_lung_samples/data_C54ctr_cb.rds")
SaveH5Seurat(c54, filename = "data/Seurat2AdataIntermediateFiles/C54_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C54_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c54, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c54.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

c55 <- readRDS("data/raw_data/Control_lung_samples/data_C55ctr_cb.rds")
SaveH5Seurat(c55, filename = "data/Seurat2AdataIntermediateFiles/C55_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C55_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c55, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c55.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

c56 <- readRDS("data/raw_data/Control_lung_samples/data_C56ctr_cb.rds")
SaveH5Seurat(c56, filename = "data/Seurat2AdataIntermediateFiles/C56_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C56_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c56, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c56.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


c57 <- readRDS("data/raw_data/Control_lung_samples/data_C57ctr_cb.rds")
SaveH5Seurat(c57, filename = "data/Seurat2AdataIntermediateFiles/C57_ctrl.h5Seurat")
Convert("data/Seurat2AdataIntermediateFiles/C57_ctrl.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = c57, slot = "counts")), 
            'data/Seurat2AdataIntermediateFiles/ctrl_c57.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)




