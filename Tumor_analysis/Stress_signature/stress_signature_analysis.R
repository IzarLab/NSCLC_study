# Author: Somnath Tagore, Ph.D. 
# Title: Cluster 21 analysis 
# Script Name: cluster_21_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs analysis on Cluster 21
#
# -------------------------
#!/usr/bin/env Rscript

# Load packages
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')

# 1) Stress_module Primary vs BM
# -----------------------------
#

#manual color
ichorcna_data <- read.csv(file="./Integrated_NSCLC_tumor_stress_module.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CX17", title = "CX17: PRIMARY vs BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, size = 0.05,
                   font.label = list(size = 42, face = "bold"), 
                   ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#000000','PRIMARY'='#000000')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))

ggsave("./Stress_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)

