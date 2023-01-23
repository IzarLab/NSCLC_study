# Author: Somnath Tagore, Ph.D. Title: Slideseq rctd annotations
# Script Name: annotations_rctd.R 
# Last Updated: 10/01/2022

library(spacexr)
library(Matrix)
library(stringr)
library(Seurat)
library(ggplot2)


datadir = "/home/ubuntu/NSCLC/slide-seq/Puck_220831_34_STK_20"

##reference <- Reference(counts, cell_types)
#reference <- readRDS(file="../reference.slideseq.finer.rds")
reference <- readRDS(file="../reference.slideseq.rds")

  pat <- "Puck_34_STK_20"

  data_pat = readRDS(paste0(datadir,"/",pat,".rds"))
#  system(paste0("rm ",datadir,"/",pat,".rds"))

  #Extract coordinates and counts information for each puck
  DefaultAssay(data_pat) = "Spatial"
  barcodes = colnames(data_pat)
  genenames = rownames(data_pat)
  coords = data_pat$image@coordinates[colnames(data_pat),]

  counts = data_pat@assays$Spatial@counts
  rownames(counts) = genenames
  colnames(counts) = barcodes

  #run RCTD using single nuclei reference, and puck count and coordinate information
  #save in RDS file
  #output summary pdf figures
 # puck <- SpatialRNA(coords, counts)
 # myRCTD <- create.RCTD(puck, reference, max_cores = 4)
#  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

 # saveRDS(myRCTD,paste0(datadir,"/",pat,"_rctd_main_finer.rds"))

  
#myRCTD<-readRDS(file=paste0(datadir,"/",pat,"_rctd_main_finer.1.rds"))
myRCTD<-readRDS(file=paste0(datadir,"/",pat,"_rctd_main.rds"))
#myRCTD<-readRDS(file=paste0(datadir,"/",pat,"_rctd_main_finer.rds"))
    
  results <- myRCTD@results

 table(results$results_df$first_type)
  #table(results$results_df$second_type)


 # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- datadir


cols <- c("B_cells"= "blue",
          "CNS" = "green","Endothelial"= "yellow","Fibroblasts" ="seagreen3",
          "Myeloid"= "darkviolet",
          "T_cells"= "orangered4","Tumor"= "plum")

  # make the plots 
pdf("all.pdf")
  # Plots the confident weights for each cell type as in full_mode (saved as 
  # 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots the weights for each cell type as in doublet_mode. (saved as 
  # 'results/cell_type_weights_doublets.pdf')
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
		       results$results_df) 
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
  # 'results/cell_type_occur.pdf')
 plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)+scale_fill_manual(values=cols)+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

  # makes a map of all cell types, (saved as 
  # 'results/all_cell_types.pdf')
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)+scale_color_manual(values= cols)


  # doublets
  #obtain a dataframe of only doublets
  doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
  # Plots all doublets in space (saved as 
  # 'results/all_doublets.pdf')
  plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names)+scale_color_manual(values= cols)

  # Plots all doublets in space for each cell type (saved as 
  # 'results/all_doublets_type.pdf')
  plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names)+scale_color_manual(values= cols)
  # a table of frequency of doublet pairs 
  doub_occur <- table(doublets$second_type, doublets$first_type) 
  # Plots a stacked bar plot of doublet ocurrences (saved as 
  # 'results/doublet_stacked_bar.pdf')

  plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names)
#}
dev.off()
