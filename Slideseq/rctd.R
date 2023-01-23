# Author: Somnath Tagore, Ph.D. Title: Slideseq rctd analysis
# Script Name: rctd.R 
# Last Updated: 10/01/2022



library(spacexr)
library(Matrix)
library(stringr)
library(Seurat)

### title: Use RCTD pipeline to assign cell type identities to SlideSeq samples
### author: Yiping Wang date: 03/29/2022

datadir = "/home/ubuntu/NSCLC/slide-seq/Puck_220831_23_PA056"

reference <- readRDS(file="reference.slideseq.rds")
pat <- "Puck_23_PA056"
  data_pat = readRDS(paste0(datadir,"/",pat,".rds"))
  system(paste0("rm ",datadir,"/",pat,".rds"))

  #Extract coordinates and counts information for each puck
  DefaultAssay(data_pat) = "Spatial"
  barcodes = colnames(data_pat)
  genenames = rownames(data_pat)
  coords = data_pat$image@coordinates[colnames(data_pat),]
  # if (pat=="puck5" || pat=="puck7_20_feature_threshold" || pat=="puck8_20_feature_threshold")
  # {
  #   coords = coords[,c("x","y")]
  # }
  # if (pat=="puck6final")
  # {
  #   #coords = data.frame(x=coords$xcoord,y=coords$ycoord)
  #   coords$x = coords$xcoord
  #   coords$y = coords$ycoord
  #   coords$xcoord = NULL
  #   coords$ycoord = NULL
  # }
  counts = data_pat@assays$Spatial@counts
  rownames(counts) = genenames
  colnames(counts) = barcodes

  #run RCTD using single nuclei reference, and puck count and coordinate information
  #save in RDS file
  #output summary pdf figures
  puck <- SpatialRNA(coords, counts)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

  saveRDS(myRCTD,paste0(datadir,"/",pat,"_rctd_main.rds"))
