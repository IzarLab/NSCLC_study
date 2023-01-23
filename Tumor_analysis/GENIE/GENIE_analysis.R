library(matrixStats)
library(copynumber)
library(infercnv)
library(stringr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
# Author: Somnath Tagore, Ph.D. 
# Title: GENIE analysis 
# Script Name: GENIE_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs Fraction of Genome Altered in GENIE cohort 
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

# 1) GENIE NSCLC All STK_mut vs All other genes
# -----------------------------
ichorcna_data <- read.csv(file="GENIE_NSCLC_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)

# Order boxes by median
group_ordered <- with(ichorcna_data,                       
                      reorder(STK11_MUT_vs_Other_genes,
                              Fraction.Genome.Altered,
                              median))
# Print order
group_ordered                                     
data_ordered <- ichorcna_data                              
# Create data with reordered group levels
data_ordered$STK11_MUT_vs_Other_genes <- factor(data_ordered$STK11_MUT_vs_Other_genes,
                                                levels = levels(group_ordered))

yplot <- ggboxplot(data_ordered, x = "STK11_MUT_vs_Other_genes", y = "Fraction.Genome.Altered", title = "GENIE (NSCLC): STK11-MUT (All samples) vs STK11-WT (Other gene mutations)",
                   color = "STK11_MUT_vs_Other_genes", fill = "STK11_MUT_vs_Other_genes", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20)) #+ scale_fill_manual(breaks = c('EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#values=c('#696969','#3CB371','#E7BB00','#4B0082')) 
yplot + theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11_MUT_vs_Other_genes)),
           y = aggregate(Fraction.Genome.Altered ~ STK11_MUT_vs_Other_genes, data_ordered, median)[ , 2],
           label = table(data_ordered$STK11_MUT_vs_Other_genes),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values = c('EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082')) + scale_fill_manual(values=c('EGFR_MUT'='#696969','STK11_MUT'='#000000','KRAS_MUT'='#E7BB00','Other_genes_MUT'='#4B0082'))
ggsave("Fraction_Genome_Altered_STK11_NSCLC_GENIE_other_genes_v3.pdf",height = 10,width = 10)

# 2) GENIE NSCLC All STK11_mut vs STK11_wt
# -----------------------------

ichorcna_data <- read.csv(file="GENIE_NSCLC_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "Fraction.Genome.Altered", title = "GENIE (NSCLC): STK-mut (All samples) vs Non-STK-mut (All samples)",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
ggsave("Fraction_Genome_Altered_STK11_NonSTK11_NSCLC_GENIE_v2.pdf",height = 15,width = 15)

# 3) GENIE LUAD All STK_mut vs All other genes
# -----------------------------

ichorcna_data <- read.csv(file="GENIE_LUAD_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Other_genes", y = "Fraction.Genome.Altered", title = "GENIE (LUAD): STK-mut (All samples) vs Non-STK-mut (Other gene mutations)",
                   color = "STK11_mut_vs_Other_genes", fill = "STK11_mut_vs_Other_genes", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
ggsave("Fraction_Genome_Altered_STK11_LUAD_GENIE_other_genes_v2.pdf",height = 15,width = 15)

# 4) GENIE LUAD All STK11_mut vs STK11_wt
# -----------------------------

ichorcna_data <- read.csv(file="GENIE_LUAD_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "Fraction.Genome.Altered", title = "GENIE (LUAD): STK-mut (All samples) vs Non-STK-mut (All samples)",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
ggsave("Fraction_Genome_Altered_STK11_NonSTK11_LUAD_GENIE_other_genes_v2.pdf",height = 15,width = 15)

# 5) GENIE LUSC All STK_mut vs All other genes
# -----------------------------
ichorcna_data <- read.csv(file="GENIE_LUSC_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Other_genes", y = "Fraction.Genome.Altered", title = "GENIE (LUSC): STK-mut (All samples) vs Non-STK-mut (Other gene mutations)",
                   color = "STK11_mut_vs_Other_genes", fill = "STK11_mut_vs_Other_genes", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
ggsave("Fraction_Genome_Altered_STK11_LUSC_GENIE_other_genes_v2.pdf",height = 15,width = 15)

# 6) GENIE LUSC All STK11_mut vs STK11_wt
# -----------------------------

ichorcna_data <- read.csv(file="GENIE_LUSC_Fraction_Genome_Altered.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "Fraction.Genome.Altered", title = "GENIE (LUSC): STK-mut (All samples) vs Non-STK-mut (All samples)",
                   color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
                   alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
yplot + stat_compare_means()+theme(text = element_text(size = 20))
ggsave("Fraction_Genome_Altered_STK11_NonSTK11_LUSC_GENIE_other_genes_v2.pdf",height = 15,width = 15)

