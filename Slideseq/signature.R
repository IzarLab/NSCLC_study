# Author: Somnath Tagore, Ph.D. Title: Slideseq signature analysis
# Script Name: signature.R 
# Last Updated: 10/01/2022



library(ggplot2)
library(Seurat)
#library(ggrastr)
library(rlist)
library(RColorBrewer)
library(grid)
library(stringr)
library(ggpubr)

pats = c("Puck_23_PA056")
output_pats = c("Puck_23_PA056")
output_pats_order = output_pats

sigs_to_test = list()
#cin70 <- list(c('GFAP', 'FN1', 'FLT1', 'TCF4', 'ZEB1', 'SPARCL1', 'MBP', 'LDB2', 'ZEB2',  'IGFBP7',  'PLP1', 'CACNA1C', 'NRG3', 'PCDH9', 'KCNIP4','LAMA2','MCAM','COL4A1', 'PPP2R2B','LTBP1'))
cin70 <- list(c('TPX2','PRC1','FOXM1 ','CDC2','C20orf24','TGIF2','MCM2','H2AFZ','TOP2A','PCNA ','UBE2C ','MELK','TRIP13 ','CNAP1','MCM7','RNASEH2A ','RAD51AP1 ','KIF20A','CDC45L','MAD2L1','ESPL1 ','CCNB2','FEN1 ','TTK','CCT5 ','RFC4','ATAD2 ','ch-TOG','NUP205 ','CDC20','CKS2 ','RRM2','ELAVL1','CCNB1','RRM1','AURKB ','EZH2','CTPS ','DKC1','OIP5 ','CDCA8','PTTG1 ','CEP55 ','H2AFX ','CMAS','BRRN1','MCM10 ','LSM4','MTB','ASF1B','ZWINT','TOPK','FLJ10036','CDCA3','ECT2','CDC6','UNG ','MTCH2','RAD21','ACTL6A','GPIandMGC13096','SFRS2','HDGF','NXT1','NEK2','DHCR7','STK6','NDUFAB1','KIAA0286','KIF4A'))

names(cin70) = c("test_sig")

sigs_to_test_cin70 = cin70

min_score = 100
max_score = -100

puckarr = list()
pat = pats

puck = readRDS(paste0(pat,".rds"))

DefaultAssay(puck) = "SCT"

rctd = readRDS(paste0(pat,"_rctd_main.rds"))
puck$rctd_cell_type = ""
puck$rctd_cell_type[match(rownames(rctd@results$results_df),colnames(puck))] = as.vector(rctd@results$results_df$first_type)
puck = subset(puck, rctd_cell_type=="Tumor")

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
SpatialColors <- colorRampPalette(colors = brewer.pal(n = 24, name = "Blues"))

breaks_arr = seq(from=min_score[[1]], to=max_score[[1]], length.out=100)
colors_arr = SpatialColors(n=100)

pdf(paste0("plot_signatures_pub_quality_",pat,"_cin70.pdf"))
theme_set(theme_bw())
for (z in 1:length(pats)) {
  puck = puckarr[[pats[z]]]
 for (z1 in 1:length(names(sigs_to_test_cin70)))
 {
    sig_to_test = names(sigs_to_test_cin70)[z1]
    plotdf = data.frame(x=puck$image@coordinates[colnames(puck),]$x, y=puck$image@coordinates[colnames(puck),]$y,
                        color = puck[[paste0(sig_to_test,"1")]][[1]])
    write.csv(plotdf,file=paste0("plotdf_",pat,"_cin70.csv"))
    plotdf$xhold = plotdf$x
    plotdf$x = plotdf$y
    plotdf$y = -(plotdf$xhold)
    print(ggplot(plotdf, aes(x=x, y=y, color = color)) + scale_color_gradientn(breaks = breaks_arr, colors = colors_arr) + geom_point(size=.1)+ ggtitle("Puck_220831_23_PA056: CIN70"))#rasterise(geom_point(size=.1)) +
  }
}
dev.off()


