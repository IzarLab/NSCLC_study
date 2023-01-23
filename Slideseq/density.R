# Author: Somnath Tagore, Ph.D. Title: Slideseq density contour map
# Script Name: density.R 
# Last Updated: 10/01/2022

library(ggplot2)
library(Seurat)
#library(ggrastr)
library(rlist)
library(RColorBrewer)
library(grid)
library(stringr)
library(ggpubr)
library(tidyverse)

pats = c("Puck_23_PA056")
output_pats = c("Puck_23_PA056")
output_pats_order = output_pats

sigs_to_test = list()

myeloid_sig <- list(c('ITGAM','CD14','CD68','P2RY12','CD83','CLEC7A','CD33','FCER1A','CEACAM8','IL7R'))

cin70 <- myeloid_sig
#cin70 <- tumor_sig

names(cin70) = c("test_sig")

sigs_to_test_cin70 = cin70

min_score = 100
max_score = -100

puckarr = list()
pat = pats

puck = readRDS(paste0(pat,".rds"))

DefaultAssay(puck) = "SCT"

rctd = readRDS(paste0(pat,"_rctd_main_finer.1.rds"))
puck$rctd_cell_type = ""

table(rctd@results$results_df$first_type)


table(rctd@results$results_df$first_type)
puck$rctd_cell_type[match(rownames(rctd@results$results_df),colnames(puck))] = as.vector(rctd@results$results_df$first_type)

puck = subset(puck, rctd_cell_type=="Myeloid")

for (j in 1:length(sigs_to_test_cin70)) {
    asig = sigs_to_test_cin70[[j]]
        puck = AddModuleScore(puck, features = list(na.omit(asig)), name = names(sigs_to_test_cin70)[j], assay = "SCT", search = T, nbin = 3)
    current_min = quantile(puck[[paste0(names(sigs_to_test_cin70)[j],"1")]][[1]],.1)
    if (current_min < min_score)
    {
      min_score = current_min
    }
    current_max = quantile(puck[[paste0(names(sigs_to_test_cin70)[j],"1")]][[1]],.9)
    if (current_max > max_score)
    {
      max_score = current_max
    }
}

puckarr = list.append(puckarr, puck)
names(puckarr) = pats


#SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
#SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 7, name = "RdBu")))
#SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 24, name = "RdBu")))
SpatialColors <- colorRampPalette(colors = brewer.pal(n = 7, name = "Greys"))

breaks_arr = seq(from=min_score[[1]], to=max_score[[1]], length.out=10)
colors_arr = SpatialColors(n=10)

library(cowplot)
#library(graphics)
library(ggExtra)
pdf(paste0("plot_signatures_pub_quality_",pat,"_myeloid_1.pdf"),height=7,width=7)
#pdf(paste0("plot_signatures_pub_quality_",pat,"_tumor_1.pdf"),height=7,width=7)
theme_set(theme_bw())
for (z in 1:length(pats)) {
  puck = puckarr[[pats[z]]]
 for (z1 in 1:length(names(sigs_to_test_cin70)))
 {
    sig_to_test = names(sigs_to_test_cin70)[z1]
    plotdf = data.frame(x=puck$image@coordinates[colnames(puck),]$x, y=puck$image@coordinates[colnames(puck),]$y,
                        color = puck[[paste0(sig_to_test,"1")]][[1]])
    write.csv(plotdf,file=paste0("plotdf_",pat,"_myeloid.csv"))
    plotdf$xhold = plotdf$x
    plotdf$x = plotdf$y
    plotdf$y = -(plotdf$xhold)
    sp<-ggplot(plotdf, aes(x=x, y=y, color = color)) + scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature")+geom_density_2d_filled(bins=5) + scale_fill_brewer()#geom_density_2d_filled(contour_var="ndensity",bins=5)#+ggdensity(plotdf,"x",fill="color")+ggdensity(plotdf,"y",fill="color"))#+rasterise(geom_point(size=.1)))
    #plotdf.1 <- plotdf

    ggbld <- ggplot_build(sp)
    gdata <- ggbld$data[[1]]
    gdata[['barcodes']] <- puck$image@coordinates[colnames(puck),]

    write.csv(gdata,file=paste0("plotdf_",pat,"_tumor_ggbuild.csv"))
    puck$image@coordinates[colnames(puck),]
    
    xplot <- ggdensity(plotdf,"x",fill="color")
    yplot <- ggdensity(plotdf,"y",fill="color")
    rotate()
    sp <- sp + rremove("legend")
    xplot <- xplot + clean_theme() + rremove("legend")
    yplot <- yplot + clean_theme() + rremove("legend")
    #print(plot_grid(xplot,NULL,sp,yplot,nccol=2,align="hv",rel_widths=c(2,1),rel_heights=c(1,2))) 
    #print(ggMarginal(sp,type="density",fill="grey",color="grey")) 
    print(ggMarginal(sp,type="histogram"))
   #sp1<-ggplot(plotdf, aes(x=x, y=y, color = color)) + scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature")+geom_density_2d_filled() + scale_fill_brewer() 
   #print(ggMarginal(sp1,type="histogram")) 
   #print(ggplot(plotdf) + geom_raster(aes(x = x, y = y, fill=color)))#+ scale_fill_gradientn(colors = SpatialColors,  na.value="white") + ggtitle("Puck_220831_23_PA056: Myeloid Signature"))# +  ggplot2::scale_shape_identity() + ggplot2::theme_classic() + ggplot2::scale_size_identity() + coord_fixed()
    #print(ggplot(plotdf, aes(x=x,y=y,color=color)) + geom_tile() + ggtitle("Puck_220831_23_PA056: Myeloid Signature")) 
	 # scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature"))#rasterise(geom_point(size=.1)) +
  }
}
dev.off()

pdf(paste0("plot_signatures_pub_quality_",pat,"_myeloid_2.pdf"),height=7,width=7)
#pdf(paste0("plot_signatures_pub_quality_",pat,"_tumor_2.pdf"),height=7,width=7)
theme_set(theme_bw())
for (z in 1:length(pats)) {
  puck = puckarr[[pats[z]]]
 for (z1 in 1:length(names(sigs_to_test_cin70)))
 {
    sig_to_test = names(sigs_to_test_cin70)[z1]
    plotdf = data.frame(x=puck$image@coordinates[colnames(puck),]$x, y=puck$image@coordinates[colnames(puck),]$y,
                        color = puck[[paste0(sig_to_test,"1")]][[1]])
 #   write.csv(plotdf,file=paste0("plotdf_",pat,"_tumor.csv"))
    plotdf$xhold = plotdf$x
    plotdf$x = plotdf$y
    plotdf$y = -(plotdf$xhold)
    sp<-ggplot(plotdf, aes(x=x, y=y, color = color)) + scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature")+geom_density_2d_filled(bins=5) + scale_fill_brewer()#geom_density_2d_filled(contour_var="ndensity",bins=5)#+ggdensity(plotdf,"x",fill="color")+ggdensity(plotdf,"y",fill="color"))#+rasterise(geom_point(size=.1)))
    xplot <- ggdensity(plotdf,"x",fill="color")
    yplot <- ggdensity(plotdf,"y",fill="color")
    rotate()
    sp <- sp + rremove("legend")
    xplot <- xplot + clean_theme() + rremove("legend")
    yplot <- yplot + clean_theme() + rremove("legend")
    #print(plot_grid(xplot,NULL,sp,yplot,nccol=2,align="hv",rel_widths=c(2,1),rel_heights=c(1,2))) 
    #print(ggMarginal(sp,type="density",fill="grey",color="grey")) 
#    print(ggMarginal(sp,type="histogram"))
   sp1<-ggplot(plotdf, aes(x=x, y=y, color = color)) + scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature")+geom_density_2d_filled() + scale_fill_brewer()
   sp1 <- sp1 + rremove("legend")
   print(ggMarginal(sp1,type="histogram"))
   #print(ggplot(plotdf) + geom_raster(aes(x = x, y = y, fill=color)))#+ scale_fill_gradientn(colors = SpatialColors,  na.value="white") + ggtitle("Puck_220831_23_PA056: Myeloid Signature"))# +  ggplot2::scale_shape_identity() + ggplot2::theme_classic() + ggplot2::scale_size_identity() + coord_fixed()
    #print(ggplot(plotdf, aes(x=x,y=y,color=color)) + geom_tile() + ggtitle("Puck_220831_23_PA056: Myeloid Signature")) 
         # scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: Myeloid Signature"))#rasterise(geom_point(size=.1)) +
  }
}
dev.off()


