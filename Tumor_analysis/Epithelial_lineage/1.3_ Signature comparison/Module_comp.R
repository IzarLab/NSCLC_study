library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)
library(celldex)
library(SingleR)



all_data <- readRDS("../Integrated_NSCLC_KRAS_STK_samples_tumor_44_samples_annotations.rds")

print(all_data)

DefaultAssay(all_data) <- "RNA"


modules <- read.csv(file = 'module_marker_list.csv')
for (i in 1:length(modules)){
    genes <- (lapply(modules[i], function(x) x[x != ""]))[[1]]
    name <- colnames(modules)[i]   
    all_data<-AddModuleScore(all_data, list(genes), name = name)
    write.csv(all_data@meta.data,"modulescores.csv")

    print("added_one")                 
}
write.csv(all_data@meta.data,"modulescores.csv")
