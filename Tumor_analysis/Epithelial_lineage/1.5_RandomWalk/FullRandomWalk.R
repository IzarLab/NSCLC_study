library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)
library(celldex)
library(SingleR)
library(Seurat)
library(dplyr)
library(pheatmap)

rds_obj <-readRDS('../TumorSampleIntegration/Full500MergedTumorCtrlIntegratedAssay.rds')
data <- read.csv("all_cell_assignments.csv")

rownames(data)<- data$Unnamed..0.1
subsampled_seurat_obj_to_map <- subset(rds_obj, cells = rownames(data))
set.seed(0) # for reproducibility
sample_fraction <- 0.1 # proportion of cells to sample per group
all_rows_to_keep <- c()
# Function to subsample cells
for (x in unique(subsampled_seurat_obj_to_map$orig.ident)){
	    subs <- subset(subsampled_seurat_obj_to_map@meta.data, orig.ident == x)
    smaller_sub <- sample(nrow(subs),nrow(subs)*sample_fraction) 
        all_rows <- rownames(subs[smaller_sub, ])
        all_rows_to_keep <- c(all_rows_to_keep, all_rows)
}


data_2<-read.csv("Cluster_21.csv")
data_2$X <- paste("Tumor_", data_2$X, sep="")
data_2$X 

all_rows_to_keep<-c(all_rows_to_keep, data_2$X )
length(all_rows_to_keep)
all_rows_to_keep<- unique(all_rows_to_keep)
subsampled_seurat_obj <- subset(subsampled_seurat_obj_to_map, cells = all_rows_to_keep)
subsampled_seurat_obj <- RunPCA(subsampled_seurat_obj,
				                                features = VariableFeatures(object = subsampled_seurat_obj))

subsampled_seurat_obj <- FindNeighbors(subsampled_seurat_obj, dims = 1:10)
subsampled_seurat_obj <- FindClusters(subsampled_seurat_obj, resolution = 0.5)
dictionary <- setNames(data$Nearest_clustering_PCA, rownames(data))
subsampled_seurat_obj@meta.data$Nearest_clustering_PCA <- dictionary[rownames(subsampled_seurat_obj@meta.data)]
graph <- subsampled_seurat_obj@graphs$integrated_snn
transition_matrix <- as.matrix(as(graph, "dgCMatrix")) # Convert the graph to a matrix
transition_matrix <- transition_matrix / Matrix::rowSums(transition_matrix)
prob_matrix <- transition_matrix
unique(subsampled_seurat_obj@meta.data$Nearest_clustering_PCA)

# Function to simulate a random walk from a given starting cell
random_walk <- function(start_cell, num_steps, prob_matrix) {
	  current_cell <- start_cell
  for (i in 1:num_steps) {
	      current_cluster <- dictionary[[current_cell]]
          if (current_cluster != 'Mixed Lineage') {
		            return(current_cluster)
	          }
	      current_cell <- sample(x = colnames(prob_matrix),
				                                    size = 1,
								                                   prob = prob_matrix[current_cell, ])
	    }
    return(dictionary[[current_cell]])
}

unassigned_cells = rownames(
			    subset(subsampled_seurat_obj@meta.data,
				          Nearest_clustering_PCA == "Mixed Lineage"))
       
       
# Load the parallel package
library(parallel)
detectCores()
# Create a cluster using the number of available cores
no_cores <- 10  # leave one core free for system processes
cl <- makeCluster(no_cores)

# Export necessary objects and functions to the cluster
clusterExport(cl, varlist = c("prob_matrix", "random_walk", "dictionary"))

timing_results <- system.time({
	  results <- parSapply(cl, unassigned_cells, function(cell) {
				           final_clusters <- replicate(100, random_walk(start_cell = cell, num_steps = 50, prob_matrix))
					       table(final_clusters) / length(final_clusters)  # Normalizing to get probabilities
					     })
})

# Stop the cluster
stopCluster(cl)
saveRDS(results, "random_walk_full_c21_inclusion.rds")

print(paste("Time taken for random walks: ", timing_results['elapsed'], " seconds"))

