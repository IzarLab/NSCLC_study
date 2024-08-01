
library(Seurat)
#library(SeuratData)
library(SeuratDisk)

rds_temp <- readRDS("/home/ubuntu/durva/durva_rt_integrated_v3.2.rds")
SaveH5Seurat(rds_temp, filename = "/home/ubuntu/durva/durva_rt_integrated_v3.2.h5Seurat")
Convert("/home/ubuntu/durva/durva_rt_integrated_v3.2.h5Seurat", dest = "h5ad")

write.table(as.matrix(GetAssayData(object = rds_temp, slot = "counts")),
            '/home/ubuntu/durva/durva_rt_integrated_v3.2.csv',
            sep = ',', row.names = T, col.names = T, quote = F)

