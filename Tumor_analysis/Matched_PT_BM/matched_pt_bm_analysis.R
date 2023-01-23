# Author: Somnath Tagore, Ph.D. 
# Title: Matched PT/BM analysis 
# Script Name: matched_pt_bm_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs analysis on matched PT/BM analysis
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

# 1) Matched PT and BM
# -----------------------------

#N254 and PA060
#manual color
ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_NSCLC_tumor_1.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(#ichorcna_data, x = "Samples", y = "CIN70", title = "CIN70 Signature: N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "Cluster_21", title = "Cluster_21 Signature: N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX1", title = "CX1 (CIN.17.archetypes) Signature: N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX2", title = "CX2 (CIN.17.archetypes) Signature: N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX3", title = "CX3 (CIN.17.archetypes) Signature: N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX4", title = "CX4 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX5", title = "CX5 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX6", title = "CX6 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX7", title = "CX7 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX8", title = "CX8 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  #ichorcna_data, x = "Samples", y = "CX9", title = "CX9 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  ichorcna_data, x = "Samples", y = "CX10", title = "CX9 (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)",
  color = "Samples", fill = "Samples", palette = "jco",
  alpha = 0.5, size = 0.05,
  font.label = list(size = 42, face = "bold"), 
  ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('N254'='#000000','PA060'='#FF0000')) +scale_fill_manual(values=c('N254'='#000000','PA060'='#FF0000'))
#ggsave("CIN70_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("Cluster_21_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX1_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX2_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX3_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX4_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX5_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX6_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX7_signature_N254_PA060_v6.pdf",height = 10,width = 10)
#ggsave("CX8_signature_N254_PA060_v6.pdf",height = 10,width = 10)
ggsave("CX9_signature_N254_PA060_v6.pdf",height = 10,width = 10)


#ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_NSCLC_tumor_1.csv")
ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_2.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
for(i in 1:17)
{
  yplot <- ggboxplot(#ichorcna_data,  x = "Samples", y = paste0('CX',i), title = paste0("CX",i," (CIN.17.archetypes) Signature\n N254 (PT) and PA060 (BM)"),
    ichorcna_data,  x = "Samples", y = paste0('CX',i), title = paste0("CX",i," (CIN.17.archetypes) Signature\n N586 (PT) and PA067 (BM)"),
    color = "Samples", fill = "Samples", palette = "jco",
    alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
  #yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('N254'='#000000','PA060'='#FF0000')) +scale_fill_manual(values=c('N254'='#000000','PA060'='#FF0000'))
  yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('N586'='#A52A2A','PA067'='#191970')) +scale_fill_manual(values=c('N586'='#A52A2A','PA067'='#191970'))
  
  ggsave(paste0("CX",i,"_signature_N586_PA067_v6.pdf"),height = 10,width = 10)
}

ichorcna_data <- read.csv(file="Integrated_matched_PT_BM_2.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
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
