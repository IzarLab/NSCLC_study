# Author: Somnath Tagore, Ph.D. 
# Title: Cluster 21 analysis 
# Script Name: cluster_21_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs analysis on Cluster 21
#
# -------------------------
#!/usr/bin/env Rscript

library(matrixStats)
library(copynumber)
library(infercnv)
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(ggExtra)
library(cowplot)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(gplots)
library(singscore)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(reshape2)
library(ggpubr)

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

# 1) Stress_module Primary vs BM
# -----------------------------
#

#manual color
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_vanhove.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_marsh.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_denisenko.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_brink.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_mito_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_Ig_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_IFN_module.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/CIN_signature/Integrated_NSCLC_tumor_CX.17.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CX17", title = "CX17: PRIMARY vs BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, size = 0.05,
                   font.label = list(size = 42, face = "bold"), 
                   ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#000000','PRIMARY'='#000000')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))

#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CIN70_signature_PRIMARY_vs_BRAIN_METS_v2.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CIN70_signature_PRIMARY_vs_BRAIN_METS_v5.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_module_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_vanhove_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_marsh_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_denisenko_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_brink_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/mito_module_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/Ig_module_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
#ggsave("~/Documents/Ben_Izar_Project/stress_signature/IFN_module_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)
ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CX17_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)

# 2) Stress_module STK11mut_vs_NonSTK11mut
# -----------------------------
#

#stress_module 
#manual color
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_vanhove.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_marsh.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_denisenko.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_brink.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_mito_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_Ig_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_IFN_module.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/CIN_signature/Integrated_NSCLC_tumor_CX.17.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
for(i in 1:17)
{
  yplot <- ggboxplot(ichorcna_data, x = "STK11mut_vs_NonSTK11mut", y = paste0('CX',i), title = paste0("CX",i,": STK11-mut vs Non-STK11-mut"),
                     color = "STK11mut_vs_NonSTK11mut", fill = "STK11mut_vs_NonSTK11mut", palette = "jco",
                     alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
  yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')) +scale_fill_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000'))
  
  #yplot  + stat_compare_means(method = "t.test")
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_module_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_vanhove_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_marsh_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_denisenko_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_brink_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/mito_module_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/Ig_module_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/IFN_module_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CX1_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/CIN_signature/CX2_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)
  ggsave(paste0("~/Documents/Ben_Izar_Project/CIN_signature/CX",i,"_STK11mut_vs_NonSTK11mut_v6.pdf"),height = 10,width = 10)
}
#stress_module Primary_BM_STK_NonSTK
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

#manual color
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_vanhove.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_marsh.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_denisenko.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_stress_brink.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_mito_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_Ig_module.csv")
#ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/stress_signature/Integrated_NSCLC_tumor_IFN_module.csv")
ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/CIN_signature/Integrated_NSCLC_tumor_CX.17.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
for(i in 1:17)
{
  yplot <- ggboxplot(ichorcna_data, x = "Primary_BM_STK_NonSTK", y = paste0('CX',i), title = paste0("CX",i,": STK11-mut (Primary) vs STK11-mut (BM) \n Non-STK11-mut (Primary) vs Non-STK11-mut (BM) "),
                     color = "Primary_BM_STK_NonSTK", fill = "Primary_BM_STK_NonSTK", palette = "jco",
                     alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
  yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=colPrimvsBMvsSTKvsNonSTKmain) +scale_fill_manual(values=colPrimvsBMvsSTKvsNonSTKmain)+rotate_x_text(90)
  
  #yplot  + stat_compare_means(method = "t.test")
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_module_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_vanhove_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_marsh_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_denisenko_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/stress_brink_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/mito_module_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/Ig_module_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  #ggsave("~/Documents/Ben_Izar_Project/stress_signature/IFN_module_Primary_BM_STK_NonSTK_v6.pdf",height = 10,width = 10)
  ggsave(paste0("~/Documents/Ben_Izar_Project/CIN_signature/CX",i,"_Primary_BM_STK_NonSTK_v6.pdf"),height = 10,width = 10)
  
}
table(ichorcna_data$STK11mut_vs_NonSTK11mut)
