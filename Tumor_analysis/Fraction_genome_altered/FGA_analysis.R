library(matrixStats)
library(copynumber)
library(infercnv)
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
# Author: Somnath Tagore, Ph.D. 
# Title: Fraction of genome altered analysis 
# Script Name: FGA_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs Fraction of Genome Altered analysis 
#
# -------------------------
#!/usr/bin/env Rscript

library(ggplot2)
library("ggpubr")
library("ggExtra")

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

# 1) FGA NSCLC Izar lab cohort (inferCNV CNA)
# -----------------------------

all_samples_table = data.frame(sampleID = character(), chrom = integer(), start.pos = integer(), end.pos = integer(), mean = double())

infercnv_pats = c('NSCLC_PA076_cnv')
#infercnv_pats = c('NSCLC_KRAS_4_cnv','NSCLC_KRAS_6_cnv','NSCLC_KRAS_7_cnv','NSCLC_KRAS_8_cnv','NSCLC_KRAS_10_cnv','NSCLC_KRAS_11_cnv','NSCLC_KRAS_12_cnv','NSCLC_KRAS_13_cnv','NSCLC_KRAS_17_cnv')
# infercnv_pats = c('NSCLC_PA001_cnv','NSCLC_PA004_cnv','NSCLC_PA005_cnv','NSCLC_PA019_cnv','NSCLC_PA025_cnv','NSCLC_PA034_cnv','NSCLC_PA042_cnv','NSCLC_PA043_cnv','NSCLC_PA048_cnv','NSCLC_PA054_cnv','NSCLC_PA056_cnv',
#                   'NSCLC_PA060_cnv','NSCLC_PA067_cnv','NSCLC_PA068_cnv','NSCLC_PA070_cnv','NSCLC_PA072_cnv','NSCLC_PA076_cnv','NSCLC_PA080_cnv','NSCLC_PA104_cnv','NSCLC_PA125_cnv','NSCLC_PA141_cnv','NSCLC_N245_cnv',
#                   'NSCLC_N561_cnv','NSCLC_N586_cnv')
#infercnv_pats = c('NSCLC_STK_14_cnv','NSCLC_STK_15_cnv','NSCLC_STK_18_cnv','NSCLC_STK_20_cnv','NSCLC_STK_21_cnv','NSCLC_STK_2_cnv','NSCLC_STK_1_cnv','NSCLC_STK_5dot1_cnv','NSCLC_STK_5dot2_cnv','NSCLC_STK_22dot2_cnv','NSCLC_STK_3_cnv')
#pdf("copynumber_plot_STK11_mut.pdf")
#pdf("copynumber_plot_Non_STK11_mut.pdf",height=15,width=15)
pdf("copynumber_plot_STK11_WT.pdf",height=15,width=15)
for (z in 1:length(infercnv_pats))
{
  #infercnv_obj@observation_grouped_cell_indices$`Epithelial cells`
  # infercnv_obj = readRDS("NSCLC_PA076_cnv.rds")
  infercnv_obj = readRDS(paste0(infercnv_pats[z],".rds"))
  #infercnv_obj = readRDS(paste0("NSCLC_STK_14_cnv.rds"))
  rowmediansarr = rowMedians(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices$`Epithelial cells`])
  rowmediansarr = rowMedians(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices$`NA`])
  rowmediansarr <- round(rowmediansarr,3)
  atable = data.frame(sampleID = rep(infercnv_pats[z], length(rowmediansarr)), 
                      chrom = substring(infercnv_obj@gene_order$chr, 4), start.pos = infercnv_obj@gene_order$start, 
                      end.pos = infercnv_obj@gene_order$stop, mean = rowmediansarr)
  all_samples_table = rbind(all_samples_table, atable)
  all_samples_table <- atable
  # all_samples_table$mean
  #write.csv(all_samples_table,file=paste0("all_samples_table_copynumber_plot_STK11_mut_",infercnv_pats[z],".csv"))
  #write.csv(all_samples_table,file=paste0("all_samples_table_copynumber_plot_Non_STK11_mut_",infercnv_pats[z],".csv"))
  write.csv(all_samples_table,file=paste0("all_samples_table_copynumber_plot_STK11_WT_",infercnv_pats[z],".csv"))
  
  all_samples_table$mean[is.nan(all_samples_table$mean)] <- 0
  
  all_samples_tablemean = table(all_samples_table$mean)
  all_cnvs = as.double(names(all_samples_tablemean))
  common_cnvs = as.double(names(all_samples_tablemean[all_samples_tablemean > 100]))
  posdists = all_cnvs - common_cnvs[2]
  posdists = posdists[posdists > 0]
  negdists = all_cnvs - common_cnvs[1]
  negdists = abs(negdists[negdists < 0])
  leastdist = min(c(posdists, negdists))
  # loss_cutoff = min(common_cnvs) - leastdist/2
  #  gain_cutoff = max(common_cnvs) + leastdist/2
  
  loss_cutoff = quantile(all_samples_table$mean, .05)
  gain_cutoff = quantile(all_samples_table$mean, .95)
  
  #gain_cutoff <- 0.83
  #loss_cutoff <- 0.23
  print(plotFreq(segments = all_samples_table, thres.gain = gain_cutoff, thres.loss = loss_cutoff, ylim = c(-5, 5), main = paste0(infercnv_pats[z], " cutoffs: ", loss_cutoff, " ", gain_cutoff)))
}
dev.off()

#read cnv files
#sample_nsclc <- read.csv(file="all_samples_table_copynumber_plot_STK11_WT_NSCLC_N586_cnv.csv")
#sample_nsclc <- read.csv(file="all_samples_table_copynumber_plot_STK11_WT_NSCLC_PA076_cnv.csv")
sample_nsclc <- read.csv(file="all_samples_table_copynumber_plot_Non_STK11_mut_NSCLC_KRAS_13_cnv.csv")
#sample_nsclc <- read.csv(file="all_samples_table_copynumber_plot_STK11_WT_NSCLC_N561_cnv.csv")
sample_nsclc <- as.data.frame(sample_nsclc)

#assign gain/loss cutoff
#gain_cutoff <- 1.051
#loss_cutoff <- 0.968

loss_cutoff = quantile(sample_nsclc$mean, .05)
loss_cutoff
gain_cutoff = quantile(sample_nsclc$mean, .95)
gain_cutoff

#sample_nsclc$chrom=='1'

#calculate chromosome-wise alterations
chr1 <- sample_nsclc[sample_nsclc$chrom=='1' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr1.altered <- sum(chr1$end.pos-chr1$start.pos)/248956422

chr2 <- sample_nsclc[sample_nsclc$chrom=='2' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr2.altered <- sum(chr2$end.pos-chr2$start.pos)/242193529

chr3 <- sample_nsclc[sample_nsclc$chrom=='3' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr3.altered <- sum(chr3$end.pos-chr3$start.pos)/198295559

chr4 <- sample_nsclc[sample_nsclc$chrom=='4' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr4.altered <- sum(chr4$end.pos-chr4$start.pos)/190214555

chr5 <- sample_nsclc[sample_nsclc$chrom=='5' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr5.altered <- sum(chr5$end.pos-chr5$start.pos)/181538259

chr6 <- sample_nsclc[sample_nsclc$chrom=='6' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr6.altered <- sum(chr6$end.pos-chr6$start.pos)/170805979

chr7 <- sample_nsclc[sample_nsclc$chrom=='7' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr7.altered <- sum(chr7$end.pos-chr7$start.pos)/159345973

chr8 <- sample_nsclc[sample_nsclc$chrom=='8' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr8.altered <- sum(chr8$end.pos-chr8$start.pos)/145138636

chr9 <- sample_nsclc[sample_nsclc$chrom=='9' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr9.altered <- sum(chr9$end.pos-chr9$start.pos)/138394717

chr10 <- sample_nsclc[sample_nsclc$chrom=='10' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr10.altered <- sum(chr10$end.pos-chr10$start.pos)/133797422

chr11 <- sample_nsclc[sample_nsclc$chrom=='11' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr11.altered <- sum(chr11$end.pos-chr11$start.pos)/135086622

chr12 <- sample_nsclc[sample_nsclc$chrom=='12' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr12.altered <- sum(chr12$end.pos-chr12$start.pos)/133275309

chr13 <- sample_nsclc[sample_nsclc$chrom=='13' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr13.altered <- sum(chr13$end.pos-chr13$start.pos)/114364328

chr14 <- sample_nsclc[sample_nsclc$chrom=='14' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr14.altered <- sum(chr14$end.pos-chr14$start.pos)/107043718

chr15 <- sample_nsclc[sample_nsclc$chrom=='15' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr15.altered <- sum(chr15$end.pos-chr15$start.pos)/101991189

chr16 <- sample_nsclc[sample_nsclc$chrom=='16' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr16.altered <- sum(chr16$end.pos-chr16$start.pos)/90338345

chr17 <- sample_nsclc[sample_nsclc$chrom=='17' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr17.altered <- sum(chr17$end.pos-chr17$start.pos)/83257441

chr18 <- sample_nsclc[sample_nsclc$chrom=='18' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr18.altered <- sum(chr18$end.pos-chr18$start.pos)/80373285

chr19 <- sample_nsclc[sample_nsclc$chrom=='19' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr19.altered <- sum(chr19$end.pos-chr19$start.pos)/58617616

chr20 <- sample_nsclc[sample_nsclc$chrom=='20' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr20.altered <- sum(chr20$end.pos-chr20$start.pos)/64444167

chr21 <- sample_nsclc[sample_nsclc$chrom=='21' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr21.altered <- sum(chr21$end.pos-chr21$start.pos)/46709983

chr22 <- sample_nsclc[sample_nsclc$chrom=='22' & (sample_nsclc$mean>gain_cutoff | sample_nsclc$mean<loss_cutoff),]
chr22.altered <- sum(chr22$end.pos-chr22$start.pos)/50818468

fraction_genome_altered <- (chr1.altered + chr2.altered + chr3.altered + chr4.altered + chr5.altered +chr6.altered+chr7.altered+chr8.altered+chr9.altered+chr10.altered+chr11.altered+chr12.altered+chr13.altered+chr14.altered+chr15.altered+chr16.altered+chr17.altered+chr18.altered
                            +chr19.altered+chr20.altered+chr21.altered+chr22.altered)/22
fraction_genome_altered


###########################################################
# 2) calculate fga for stk_wt
# -----------------------------
#STK11_WT
fraction_genome_altered_PA001 <- fraction_genome_altered
fraction_genome_altered_PA004 <- fraction_genome_altered
fraction_genome_altered_PA005 <- fraction_genome_altered
fraction_genome_altered_PA019 <- fraction_genome_altered
fraction_genome_altered_PA025 <- fraction_genome_altered
fraction_genome_altered_PA034 <- fraction_genome_altered
fraction_genome_altered_PA042 <- fraction_genome_altered
fraction_genome_altered_PA043 <- fraction_genome_altered
fraction_genome_altered_PA048 <- fraction_genome_altered
fraction_genome_altered_PA054 <- fraction_genome_altered
fraction_genome_altered_PA056 <- fraction_genome_altered
fraction_genome_altered_PA060 <- fraction_genome_altered
fraction_genome_altered_PA067 <- fraction_genome_altered
fraction_genome_altered_PA068 <- fraction_genome_altered
fraction_genome_altered_PA070 <- fraction_genome_altered
fraction_genome_altered_PA072 <- fraction_genome_altered
fraction_genome_altered_PA080 <- fraction_genome_altered
fraction_genome_altered_PA104 <- fraction_genome_altered
fraction_genome_altered_PA125 <- fraction_genome_altered
fraction_genome_altered_PA141 <- fraction_genome_altered
fraction_genome_altered_N254 <- fraction_genome_altered
fraction_genome_altered_N561 <- fraction_genome_altered
fraction_genome_altered_N586 <- fraction_genome_altered


#Non_STK11_mut
fraction_genome_altered_KRAS_4 <- fraction_genome_altered
fraction_genome_altered_KRAS_6 <- fraction_genome_altered
fraction_genome_altered_KRAS_7 <- fraction_genome_altered
fraction_genome_altered_KRAS_8 <- fraction_genome_altered
fraction_genome_altered_KRAS_10 <- fraction_genome_altered
fraction_genome_altered_KRAS_11 <- fraction_genome_altered
fraction_genome_altered_KRAS_12 <- fraction_genome_altered
fraction_genome_altered_KRAS_13 <- fraction_genome_altered
fraction_genome_altered_KRAS_17 <- fraction_genome_altered

#STK11_mut
fraction_genome_altered_STK_14 <- fraction_genome_altered
fraction_genome_altered_STK_15 <- fraction_genome_altered
fraction_genome_altered_STK_18 <- fraction_genome_altered
fraction_genome_altered_STK_20 <- fraction_genome_altered
fraction_genome_altered_STK_21 <- fraction_genome_altered
fraction_genome_altered_STK_2 <- fraction_genome_altered
fraction_genome_altered_STK_1 <- fraction_genome_altered
fraction_genome_altered_STK_5dot1 <- fraction_genome_altered
fraction_genome_altered_STK_5dot2 <- fraction_genome_altered
fraction_genome_altered_STK_22dot2 <- fraction_genome_altered
fraction_genome_altered_STK_3 <- fraction_genome_altered

# chr3 <- all_samples_table_copynumber_plot_STK11_mut_NSCLC_STK_15[all_samples_table_copynumber_plot_STK11_mut_NSCLC_STK_15$chrom=='3' & (all_samples_table_copynumber_plot_STK11_mut_NSCLC_STK_15$mean>1.051 | all_samples_table_copynumber_plot_STK11_mut_NSCLC_STK_15$mean<0.96),]
# chr3.altered <- sum(chr3$end.pos-chr3$start.pos)/242193529

# create a dataset
data <- data.frame(
  Sample_groups=c( rep("STK11_mut",11),c(rep("Non_STK11_mut",9)),rep("STK11_WT",23)),
  value=c(c(fraction_genome_altered_STK_14,
            fraction_genome_altered_STK_15,
            fraction_genome_altered_STK_18,
            fraction_genome_altered_STK_20,
            fraction_genome_altered_STK_21,
            fraction_genome_altered_STK_2,
            fraction_genome_altered_STK_1,
            fraction_genome_altered_STK_5dot1,
            fraction_genome_altered_STK_5dot2,
            fraction_genome_altered_STK_22dot2,
            fraction_genome_altered_STK_3 ),c(fraction_genome_altered_KRAS_4,
                                              fraction_genome_altered_KRAS_6,
                                              fraction_genome_altered_KRAS_7 ,
                                              fraction_genome_altered_KRAS_8 ,
                                              fraction_genome_altered_KRAS_10 ,
                                              fraction_genome_altered_KRAS_11 ,
                                              fraction_genome_altered_KRAS_12 ,
                                              fraction_genome_altered_KRAS_13 ,
                                              fraction_genome_altered_KRAS_17),c(fraction_genome_altered_PA001,
                                                                                 fraction_genome_altered_PA004,
                                                                                 fraction_genome_altered_PA005 ,
                                                                                 fraction_genome_altered_PA019 ,
                                                                                 fraction_genome_altered_PA025 ,
                                                                                 fraction_genome_altered_PA034 ,
                                                                                 fraction_genome_altered_PA042 ,
                                                                                 fraction_genome_altered_PA043 ,
                                                                                 fraction_genome_altered_PA048 ,
                                                                                 fraction_genome_altered_PA054 ,
                                                                                 fraction_genome_altered_PA056 ,
                                                                                 fraction_genome_altered_PA060 ,
                                                                                 fraction_genome_altered_PA067,
                                                                                 fraction_genome_altered_PA068 ,
                                                                                 fraction_genome_altered_PA070 ,
                                                                                 fraction_genome_altered_PA072 ,
                                                                                 fraction_genome_altered_PA080 ,
                                                                                 fraction_genome_altered_PA104 ,
                                                                                 fraction_genome_altered_PA125 ,
                                                                                 fraction_genome_altered_PA141,
                                                                                 fraction_genome_altered_N254 ,
                                                                                 fraction_genome_altered_N561 ,
                                                                                 fraction_genome_altered_N586))
)

# create a dataset
data <- data.frame(
  Sample_groups=c( rep("STK11_mut",11),c(rep("Rest",32))),
  value=c(c(fraction_genome_altered_STK_14,
            fraction_genome_altered_STK_15,
            fraction_genome_altered_STK_18,
            fraction_genome_altered_STK_20,
            fraction_genome_altered_STK_21,
            fraction_genome_altered_STK_2,
            fraction_genome_altered_STK_1,
            fraction_genome_altered_STK_5dot1,
            fraction_genome_altered_STK_5dot2,
            fraction_genome_altered_STK_22dot2,
            fraction_genome_altered_STK_3 ),c(fraction_genome_altered_KRAS_4,
                                              fraction_genome_altered_KRAS_6,
                                              fraction_genome_altered_KRAS_7 ,
                                              fraction_genome_altered_KRAS_8 ,
                                              fraction_genome_altered_KRAS_10 ,
                                              fraction_genome_altered_KRAS_11 ,
                                              fraction_genome_altered_KRAS_12 ,
                                              fraction_genome_altered_KRAS_13 ,
                                              fraction_genome_altered_KRAS_17,fraction_genome_altered_PA001,
                                              fraction_genome_altered_PA004,
                                              fraction_genome_altered_PA005 ,
                                              fraction_genome_altered_PA019 ,
                                              fraction_genome_altered_PA025 ,
                                              fraction_genome_altered_PA034 ,
                                              fraction_genome_altered_PA042 ,
                                              fraction_genome_altered_PA043 ,
                                              fraction_genome_altered_PA048 ,
                                              fraction_genome_altered_PA054 ,
                                              fraction_genome_altered_PA056 ,
                                              fraction_genome_altered_PA060 ,
                                              fraction_genome_altered_PA067,
                                              fraction_genome_altered_PA068 ,
                                              fraction_genome_altered_PA070 ,
                                              fraction_genome_altered_PA072 ,
                                              fraction_genome_altered_PA080 ,
                                              fraction_genome_altered_PA104 ,
                                              fraction_genome_altered_PA125 ,
                                              fraction_genome_altered_PA141,
                                              fraction_genome_altered_N254 ,
                                              fraction_genome_altered_N561 ,
                                              fraction_genome_altered_N586))
)

# Plot
#pdf(file='Frac_of_gen_alt_STK11_mut_vs_Non_STK11_mut_STK11_WT.pdf',height=10,width=10)
pdf(file='Frac_of_gen_alt_STK11_mut_vs_Rest.pdf',height=10,width=10)
data %>%
  ggplot( aes(x=Sample_groups, y=value, fill=Sample_groups)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  #ggtitle("Fraction of genome altered: STK11_mut vs Non_STK11_mut, STK11_WT") +
  ggtitle("Fraction of genome altered: STK11_mut vs Rest") +
  xlab("Sample groups")+ylab("Average Fraction of genome altered")+theme_bw()
#ggsave('Frac_of_gen_alt_STK11_mut_vs_Non_STK11_mut_STK11_WT.pdf')
dev.off()

# 3) ichorCNA fraction of genome altered
# -----------------------------

sample_nsclc <- read.csv(file="N561_ichorCNA.csv")
sample_nsclc <- as.data.frame(sample_nsclc)

#gain_cutoff <- 1.051
#loss_cutoff <- 0.968

loss_cutoff = quantile(sample_nsclc$median, .05)
loss_cutoff
gain_cutoff = quantile(sample_nsclc$median, .95)
gain_cutoff

chr1 <- sample_nsclc[sample_nsclc$chr=='1' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr1.altered <- sum(chr1$end-chr1$start)/248956422

chr2 <- sample_nsclc[sample_nsclc$chr=='2' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr2.altered <- sum(chr2$end-chr2$start)/242193529

chr3 <- sample_nsclc[sample_nsclc$chr=='3' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr3.altered <- sum(chr3$end-chr3$start)/198295559

chr4 <- sample_nsclc[sample_nsclc$chr=='4' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr4.altered <- sum(chr4$end-chr4$start)/190214555

chr5 <- sample_nsclc[sample_nsclc$chr=='5' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr5.altered <- sum(chr5$end-chr5$start)/181538259

chr6 <- sample_nsclc[sample_nsclc$chr=='6' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr6.altered <- sum(chr6$end-chr6$start)/170805979

chr7 <- sample_nsclc[sample_nsclc$chr=='7' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr7.altered <- sum(chr7$end-chr7$start)/159345973

chr8 <- sample_nsclc[sample_nsclc$chr=='8' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr8.altered <- sum(chr8$end-chr8$start)/145138636

chr9 <- sample_nsclc[sample_nsclc$chr=='9' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr9.altered <- sum(chr9$end-chr9$start)/138394717

chr10 <- sample_nsclc[sample_nsclc$chr=='10' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr10.altered <- sum(chr10$end-chr10$start)/133797422

chr11 <- sample_nsclc[sample_nsclc$chr=='11' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr11.altered <- sum(chr11$end-chr11$start)/135086622

chr12 <- sample_nsclc[sample_nsclc$chr=='12' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr12.altered <- sum(chr12$end-chr12$start)/133275309

chr13 <- sample_nsclc[sample_nsclc$chr=='13' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr13.altered <- sum(chr13$end-chr13$start)/114364328

chr14 <- sample_nsclc[sample_nsclc$chr=='14' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr14.altered <- sum(chr14$end-chr14$start)/107043718

chr15 <- sample_nsclc[sample_nsclc$chr=='15' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr15.altered <- sum(chr15$end-chr15$start)/101991189

chr16 <- sample_nsclc[sample_nsclc$chr=='16' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr16.altered <- sum(chr16$end-chr16$start)/90338345

chr17 <- sample_nsclc[sample_nsclc$chr=='17' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr17.altered <- sum(chr17$end-chr17$start)/83257441

chr18 <- sample_nsclc[sample_nsclc$chr=='18' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr18.altered <- sum(chr18$end-chr18$start)/80373285

chr19 <- sample_nsclc[sample_nsclc$chr=='19' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr19.altered <- sum(chr19$end-chr19$start)/58617616

chr20 <- sample_nsclc[sample_nsclc$chr=='20' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr20.altered <- sum(chr20$end-chr20$start)/64444167

chr21 <- sample_nsclc[sample_nsclc$chr=='21' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr21.altered <- sum(chr21$end-chr21$start)/46709983

chr22 <- sample_nsclc[sample_nsclc$chr=='22' & (sample_nsclc$median>gain_cutoff | sample_nsclc$median<loss_cutoff),]
chr22.altered <- sum(chr22$end-chr22$start)/50818468

fraction_genome_altered <- (chr1.altered + chr2.altered + chr3.altered + chr4.altered + chr5.altered +chr6.altered+chr7.altered+chr8.altered+chr9.altered+chr10.altered+chr11.altered+chr12.altered+chr13.altered+chr14.altered+chr15.altered+chr16.altered+chr17.altered+chr18.altered
                            +chr19.altered+chr20.altered+chr21.altered+chr22.altered)/22

fraction_genome_altered

# 4) ichorCNA fraction of genome altered
# -----------------------------

#compare inferCNV and ichorCNA

#
# 4.1) Primary vs BM
# -----------------------------


infercnv_vs_ichorcna <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_vs_ichorcna.csv")
infercnv_vs_ichorcna <- as.data.frame(infercnv_vs_ichorcna)

#infercnv_vs_ichorcna_prim_bm <-  infercnv_vs_ichorcna[infercnv_vs_ichorcna$Primary_vs_BM,]

library("ggpubr")
library("ggExtra")

infercnv_vs_ichorcna <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_vs_ichorcna.csv")
infercnv_vs_ichorcna <- as.data.frame(infercnv_vs_ichorcna)

pdf("~/Documents/Ben_Izar_Project/fraction_altered_genomes/fract_of_genom_alt_infercnv_ichorcna_stkmut_vs_nonstkmut_v6.pdf",height = 10,width = 15)
# Scatter plot colored by groups ("Species")
sp <- ggscatter(infercnv_vs_ichorcna, x = "inferCNV", y = "ichorCNA",add = "reg.line",
                #cor.coef = TRUE, title = "Fraction of genome altered: Primary vs BM (inferCNV vs ichorCNA)",
                cor.coef = TRUE, title = "Fraction of genome altered: STK11.mut_vs_Non.STK11.mut (inferCNV vs ichorCNA)",
                cor.method = "pearson",
                color = "STK11.mut_vs_Non.STK11.mut", palette = "jco",label = "Sample.Short.Name", repel = TRUE, 
                size = 3, alpha = 0.6, ggtheme = theme_bw())+ scale_color_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')) +scale_fill_manual(values=c('Non_STK11_mut'='#3CB371', 'STK11_mut'='#000000'))

#+stat_cor(aes(color = Primary_vs_BM), label.x = 3)              
# Marginal boxplot of x (top panel) and y (right panel)
xplot <- ggboxplot(infercnv_vs_ichorcna, x = "STK11.mut_vs_Non.STK11.mut", y = "inferCNV", 
                   color = "STK11.mut_vs_Non.STK11.mut", fill = "STK11.mut_vs_Non.STK11.mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())+ scale_color_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')) +scale_fill_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000'))+coord_flip()
yplot <- ggboxplot(infercnv_vs_ichorcna, x = "STK11.mut_vs_Non.STK11.mut", y = "ichorCNA",
                   color = "STK11.mut_vs_Non.STK11.mut", fill = "STK11.mut_vs_Non.STK11.mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())+ scale_color_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')) +scale_fill_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000'))


# Cleaning the plots
sp <- sp + rremove("legend")
#yplot <- yplot + clean_theme() + rremove("legend")
#xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv",# axis = "hv",
          rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()


pdf("~/Documents/Ben_Izar_Project/fraction_altered_genomes/fract_of_genom_alt_infercnv_ichorcna_primary_vs_bm_v6.pdf",height = 10,width = 15)
# Scatter plot colored by groups ("Species")
sp <- ggscatter(infercnv_vs_ichorcna, x = "inferCNV", y = "ichorCNA",add = "reg.line",
                cor.coef = TRUE, title = "Fraction of genome altered: Primary vs BM (inferCNV vs ichorCNA)",
                #cor.coef = TRUE, title = "Fraction of genome altered: STK11.mut_vs_Non.STK11.mut (inferCNV vs ichorCNA)",
                cor.method = "pearson",
                color = "Primary_vs_BM", palette = "jco",label = "Sample.Short.Name", repel = TRUE, 
                size = 3, alpha = 0.6, ggtheme = theme_bw())+ scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))

#+stat_cor(aes(color = Primary_vs_BM), label.x = 3)              
# Marginal boxplot of x (top panel) and y (right panel)
xplot <- ggboxplot(infercnv_vs_ichorcna, x = "Primary_vs_BM", y = "inferCNV", 
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())+ scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))+coord_flip()
yplot <- ggboxplot(infercnv_vs_ichorcna, x = "Primary_vs_BM", y = "ichorCNA",
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())+ scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))


# Cleaning the plots
sp <- sp + rremove("legend")
#yplot <- yplot + clean_theme() + rremove("legend")
#xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv",# axis = "hv",
          rel_widths = c(2, 1), rel_heights = c(1, 2))
dev.off()


#
# 4.2) STK_mut_Primary vs STK_mut_BM
# -----------------------------


#ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_prim_vs_STKmut_BM.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_prim_vs_STKmut_BM.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "Primary_vs_BM", y = "ichorCNA", title = "STK-mut (Primary) vs STK-mut (Brain Mets)",
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 24))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_prim_vs_STKmut_BM.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_prim_vs_STKmut_BM_v2.pdf",height = 10,width = 10)

#
# 4.3) Non_STK_mut_Primary vs Non_STK_mut_BM
# -----------------------------


ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_Non_STKmut_prim_vs_Non_STKmut_BM.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "Primary_vs_BM", y = "ichorCNA", title = "Non-STK-mut (Primary) vs Non-STK-mut (Brain Mets)",
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_Non_STKmut_prim_vs_Non_STKmut_BM.pdf",height = 10,width = 10)

#
# 4.4) All STK_mut vs All Non_STK_mut
# -----------------------------


ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_vs_Non_STKmut.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11.mut_vs_Non.STK11.mut", y = "ichorCNA", title = "STK-mut (All samples) vs Non-STK-mut (All samples)",
                   color = "STK11.mut_vs_Non.STK11.mut", fill = "STK11.mut_vs_Non.STK11.mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_STKmut_vs_Non_STKmut.pdf",height = 10,width = 10)


#
# 4.5) All Primary vs All Brain Mets
# -----------------------------


#ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(Primary_vs_BM,
                              ichorCNA,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$Primary_vs_BM <- factor(data_ordered$Primary_vs_BM,
                                     levels = levels(group_ordered))
table(data_ordered$Primary_vs_BM)

yplot <- ggboxplot(data_ordered, x='Primary_vs_BM', y = 'ichorCNA', title = "ichorCNA: PRIMARY vs BRAIN_METS",
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 24))
yplot + stat_compare_means()+theme(text = element_text(size = 24)) 
yplot + theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$Primary_vs_BM)),
           y = aggregate(ichorCNA ~ Primary_vs_BM, data_ordered, median)[ , 2],
           label = table(data_ordered$Primary_vs_BM),
           col = "black",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values=c('PRIMARY'='#0000FF','BRAIN_METS'='#FF0000')) +scale_fill_manual(values=c('PRIMARY'='#0000FF','BRAIN_METS'='#FF0000'))
#ggplot(tumor_nontumor, aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9',
#'MP10','MP11','MP12','MP13','MP14','MP15')), y = Number, fill = Type, group=1,label=Number)) + 
#yplot  + stat_compare_means(method = "t.test")+theme(text = element_text(size = 24))
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM_v2.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM_v3.pdf",height = 5,width = 5)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.pdf",height = 10,width = 10)

#All 
# 4.6) STK11-MUT vs STK11-WT
# -----------------------------


#ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_stk11mut_vs_stk11wt.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(STK11.MUT_vs_STK11.WT,
                              ichorCNA,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$STK11.MUT_vs_STK11.WT <- factor(data_ordered$STK11.MUT_vs_STK11.WT,
                                             levels = levels(group_ordered))
table(data_ordered$STK11.MUT_vs_STK11.WT)

yplot <- ggboxplot(data_ordered, x='STK11.MUT_vs_STK11.WT', y = 'ichorCNA', title = "ichorCNA: STK11-MUT vs STK11-WT",
                   color = "STK11.MUT_vs_STK11.WT", fill = "STK11.MUT_vs_STK11.WT", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 24))
yplot + stat_compare_means()+theme(text = element_text(size = 24)) 
yplot + theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11.MUT_vs_STK11.WT)),
           y = aggregate(ichorCNA ~ STK11.MUT_vs_STK11.WT, data_ordered, median)[ , 2],
           label = table(data_ordered$STK11.MUT_vs_STK11.WT),
           col = "black",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values=c('STK11-WT'='#3CB371', 'STK11-MUT'='#000000')) +scale_fill_manual(values=c('STK11-WT'='#3CB371', 'STK11-MUT'='#000000'))
#ggplot(tumor_nontumor, aes(x = factor(Sample, level = c('MP1','MP2','MP3','MP4','MP5','MP6','MP7','MP8','MP9',
#'MP10','MP11','MP12','MP13','MP14','MP15')), y = Number, fill = Type, group=1,label=Number)) + 
#yplot  + stat_compare_means(method = "t.test")+theme(text = element_text(size = 24))
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM_v2.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM_v3.pdf",height = 5,width = 5)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_stk11mut_vs_stk11wt_v3.pdf",height = 5,width = 5)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.pdf",height = 10,width = 10)

#
# 4.7) All Primary vs All Brain Mets (infercnv)
# -----------------------------


#ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_vs_ichorcna.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "Primary_vs_BM", y = "inferCNV", title = "Primary (All samples) vs Brain-Mets (All samples) \n inferCNV comparison",
                   color = "Primary_vs_BM", fill = "Primary_vs_BM", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 24))
#yplot  + stat_compare_means(method = "t.test")+theme(text = element_text(size = 24))
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_primary_vs_BM_v2.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.pdf",height = 10,width = 10)

#
# 4.8) All STK11-mut_vs_Non-STK11-mut (infercnv)
# -----------------------------

#
#ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_vs_ichorcna.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11.mut_vs_Non.STK11.mut", y = "inferCNV", title = "STK11-mut (All samples) vs Non-STK11-mut (All samples) \n inferCNV comparison",
                   color = "STK11.mut_vs_Non.STK11.mut", fill = "STK11.mut_vs_Non.STK11.mut", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
#yplot  + stat_compare_means(method = "t.test")+theme(text = element_text(size = 24))
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorCNA_altered/infercnv_STK11-mut_vs_Non-STK11-mut_v2.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/ichorcna_primary_vs_BM.pdf",height = 10,width = 10)

# 4.9) TCGA LUAD All STK_mut vs All Non_STK_mut
# -----------------------------


#
ichorcna_data <- read.csv(file="~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "Fraction.Genome.Altered", title = "TCGA (LUAD): STK-mut (All samples) vs Non-STK-mut (All samples)",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/som_temp_10_28_21_onwards_mcbookpro/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA.pdf",height = 10,width = 10)


# 4.10) TCGA LUAD All STK_mut vs All other genes
# -----------------------------

#
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_other_genes.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(STK11_MUT_vs_Other_genes,
                              Fraction.Genome.Altered,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$STK11_MUT_vs_Other_genes <- factor(data_ordered$STK11_MUT_vs_Other_genes,
                                                levels = levels(group_ordered))
table(data_ordered$STK11_MUT_vs_Other_genes)

yplot <- ggboxplot(ichorcna_data, x = "STK11_MUT_vs_Other_genes", y = "Fraction.Genome.Altered", title = "TCGA (LUAD): STK11-MUT (All samples) vs STK11-WT (Other gene mutations)",
                   color = "STK11_MUT_vs_Other_genes", 
                   fill = "STK11_MUT_vs_Other_genes", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw()) 
#yplot + stat_compare_means()+theme(text = element_text(size = 24))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
                                                                                       values=c('#00CED1','#696969','#000000','#E7BB00','#4B0082')) 
yplot + theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11_MUT_vs_Other_genes)),
           y = aggregate(Fraction.Genome.Altered ~ STK11_MUT_vs_Other_genes, data_ordered, median)[ , 2],
           label = table(data_ordered$STK11_MUT_vs_Other_genes),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('ALK_MUT'='#00CED1','EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082')) + scale_fill_manual(values=c('ALK_MUT'='#00CED1','EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082'))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_other_genes.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_other_genes_v2.pdf",height = 15,width =18)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_other_genes_v3.pdf",height = 10,width =10)

# 4.11) TCGA LUAD All STK_mut vs KRAS_mut
# -----------------------------


#
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_other_genes.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "KRAS_mut", y = "Fraction.Genome.Altered", title = "TCGA (LUAD): STK-mut (All samples) vs Non-STK-mut (KRAS_mut)",
                   color = "KRAS_mut", fill = "KRAS_mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_KRAS_mut.pdf",height = 10,width = 10)

# 4.12) TCGA LUAD All STK_mut vs Non-STK-mut : TMB
# -----------------------------


ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_TMB.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
table(ichorcna_data$STK11_MUT_vs_STK11_WT=='STK11_mut')
colnames(ichorcna_data)

group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(STK11_MUT_vs_STK11_WT,
                              TMB,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$STK11_MUT_vs_STK11_WT <- factor(data_ordered$STK11_MUT_vs_STK11_WT,
                                             levels = levels(group_ordered))
table(data_ordered$STK11_MUT_vs_STK11_WT)

yplot <- ggboxplot(ichorcna_data, x = "STK11_MUT_vs_STK11_WT", y = "TMB", title = "Tumor mutation burden - TCGA (LUAD)\n STK11-MUT vs STK11-WT",
                   color = "STK11_MUT_vs_STK11_WT", 
                   fill = "STK11_MUT_vs_STK11_WT", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw()) 
#yplot + stat_compare_means()+theme(text = element_text(size = 24))
yplot + stat_compare_means()+theme(text = element_text(size = 20))# + scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                    values=c('#00CED1','#696969','#000000','#E7BB00','#4B0082')) 
yplot + theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11_MUT_vs_STK11_WT)),
           y = aggregate(TMB ~ STK11_MUT_vs_STK11_WT, data_ordered, median)[ , 2],
           label = table(data_ordered$STK11_MUT_vs_STK11_WT),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000')) + scale_fill_manual(values=c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000'))

# yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "TMB", title = "Tumor mutation burden - TCGA (LUAD)\n STK-mut (All samples) vs Non-STK-mut (KRAS_mut)",
#                    color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
#                    alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
# yplot + stat_compare_means()+theme(text = element_text(size = 20))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut_v2.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut_v3.pdf",height = 7,width = 7)

#manual color
c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_TMB.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
table(ichorcna_data$STK11_MUT_vs_STK11_WT)#=='STK11_mut')
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_MUT_vs_STK11_WT", y = "TMB", title = "Tumor mutation burden - TCGA (LUAD)\n STK-mut (All samples) vs Non-STK-mut (KRAS_mut)",
                   color = "STK11_MUT_vs_STK11_WT", fill = "STK11_MUT_vs_STK11_WT", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000')) +scale_fill_manual(values=c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000'))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut_v2.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut_v5.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUAD_Non_STK11_mut_v6.pdf",height = 10,width = 10)


# 4.13) TCGA LUSC All STK_mut vs All other genes
# -----------------------------

ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUSC_TCGA_other_genes.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Other_genes", y = "Fraction.Genome.Altered", title = "TCGA (LUSC): STK-mut (All samples) vs Non-STK-mut (Other gene mutations)",
                   color = "STK11_mut_vs_Other_genes", fill = "STK11_mut_vs_Other_genes", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
#yplot  + stat_compare_means(method = "t.test")
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUSC_TCGA_other_genes.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUSC_TCGA_other_genes_v2.pdf",height = 15,width = 15)


# 4.14) TCGA LUAD All STK_mut vs Non-STK-mut : TMB
# -----------------------------

ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUSC_TCGA_MB.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "TMB", title = "Tumor mutation burden - TCGA (LUSC): STK-mut (All samples) vs Non-STK-mut (KRAS_mut)",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("~/Documents/Ben_Izar_Project/fraction_altered_genomes/TMB_STK11_LUSC_Non_STK11_mut.pdf",height = 10,width = 10)


