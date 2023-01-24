# Author: Somnath Tagore, Ph.D. 
# Title: Matched PT/BM analysis 
# Script Name: matched_pt_bm_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs analysis on matched PT/BM analysis
#
# -------------------------
#!/usr/bin/env Rscript

library(stringr)
library(tidyverse)
library(viridis)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')

# 1) Matched PT and BM
# -----------------------------

#N254 and PA060
#manual color
# Read file
ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_NSCLC_tumor_1.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# box plot
yplot <- ggboxplot(ichorcna_data, x = "Samples", y = "CIN70", title = "CIN70 Signature: N254 (PT) and PA060 (BM)",
  color = "Samples", fill = "Samples", palette = "jco",
  alpha = 0.5, size = 0.05,
  font.label = list(size = 42, face = "bold"), 
  ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('N254'='#000000','PA060'='#FF0000')) +scale_fill_manual(values=c('N254'='#000000','PA060'='#FF0000'))

#586 (PT) and PA067 (BM)
#Read file
ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_2.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

#box plot
yplot <- ggboxplot(ichorcna_data, x = "Samples", y = "CIN70", title = "CIN70 Signature: N586 (PT) and PA067 (BM)",
                   #ichorcna_data, x = "Samples", y = "Cluster_21", title = "Cluster_21 Signature: N586 (PT) and PA067 (BM)",
                   color = "Samples", fill = "Samples", palette = "jco",
                   alpha = 0.5, size = 0.05,
                   font.label = list(size = 42, face = "bold"), 
                   ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('N586'='#A52A2A','PA067'='#191970')) +scale_fill_manual(values=c('N586'='#A52A2A','PA067'='#191970'))
ggsave("CIN70_signature_N586_PA067_v6.pdf",height = 10,width = 10)
#ggsave("Cluster_21_signature_N586_PA067_v6.pdf",height = 10,width = 10)
