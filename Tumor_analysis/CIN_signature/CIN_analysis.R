
# Author: Somnath Tagore, Ph.D. 
# Title: CIN analysis 
# Script Name: CIN_analysis.R
# Last Updated: 06/24/2022

#Instructions
# The following script performs Chromosomal instability analysis 
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
library(ggplot2)
library("ggpubr")
library("ggExtra")

# Color palette
colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')

# 1) CIN25 signature PRIMARY_vs_BRAIN_METS
# -----------------------------

library(ggplot2)
ichorcna_data <- read.csv(file="CIN25.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CIN25_Signature", title = "CIN25 Signature: PRIMARY vs BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("CIN_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)


# 2) CIN25 signature STK11mut_vs_NonSTK11mut
# -----------------------------

ichorcna_data <- read.csv(file="CIN25.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "STK11mut_vs_NonSTK11mut", y = "CIN25_Signature", title = "CIN25 Signature: STK11-mut vs Non-STK11-mut",
                   color = "STK11mut_vs_NonSTK11mut", fill = "STK11mut_vs_NonSTK11mut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("CIN_signature_STK11mut_vs_NonSTK11mut.pdf",height = 10,width = 10)

colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colPrimvsBMmainTweakBM <- c('BRAIN_METS'='#FFF0F0','PRIMARY'='#0000FF')
colPrimvsBMmainTweakPRIMARY<- c('BRAIN_METS'='#FF0000','PRIMARY'='#EBF7FF')
colSTKvsNonSTKmain <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
colSTKvsNonSTKTweakSTK <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#F1F1F1')
colSTKvsNonSTKTweakNonSTK <- c('Non-STK11-mut'='#EBFFE3', 'STK11-mut'='#000000')

# 3) CIN70 signature PRIMARY_vs_BRAIN_METS
# -----------------------------

ichorcna_data <- read.csv(file="CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CIN70_Signature", title = "CIN70 Signature: PRIMARY vs BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, 
                   font.label = list(size = 42, face = "bold"), 
                   ggtheme = theme_bw())+scale_fill_manual(breaks = ichorcna_data$PRIMARY_vs_BRAIN_METS, values = c('BRAIN_METS'='black','PRIMARY'='red'))
yplot + stat_compare_means()+theme(text = element_text(size = 20))

ggplot(ichorcna_data, aes(x = PRIMARY_vs_BRAIN_METS, y = CIN70_Signature, fill = PRIMARY_vs_BRAIN_METS)) +  # Manually specified filling color
  geom_boxplot(width=0.5,lwd=1.5) +scale_fill_manual(breaks = ichorcna_data$PRIMARY_vs_BRAIN_METS, values = c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))

#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)
#ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS_v2.pdf",height = 10,width = 10)
ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS_v3.pdf",height = 10,width = 10)

#manual color
ichorcna_data <- read.csv(file="CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CIN70_Signature", title = "CIN70 Signature: PRIMARY vs BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, size = 0.05,
                   font.label = list(size = 42, face = "bold"), 
                   ggtheme = theme_bw())
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))
yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('BRAIN_METS'='#000000','PRIMARY'='#000000')) +scale_fill_manual(values=c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF'))

#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)
#ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS_v2.pdf",height = 10,width = 10)
#ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS_v5.pdf",height = 10,width = 10)
ggsave("CIN70_signature_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)


# 4) CIN archetypes signature PRIMARY_vs_BRAIN_METS
# -----------------------------

ichorcna_data <- read.csv(file="Integrated_NSCLC_tumor_CX.17.csv")
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
ggsave("CX17_PRIMARY_vs_BRAIN_METS_v6.pdf",height = 10,width = 10)

# 5) CIN70 signature STK11mut_vs_NonSTK11mut
# -----------------------------

ichorcna_data <- read.csv(file="CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(STK11_MUT_vs_STK11_WT,
                              desc(CIN70_Signature),
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$STK11_MUT_vs_STK11_WT <- factor(data_ordered$STK11_MUT_vs_STK11_WT,
                                             levels = levels(group_ordered))
table(data_ordered$STK11_MUT_vs_STK11_WT)

yplot <- ggboxplot(data_ordered, x = "STK11_MUT_vs_STK11_WT", y = "CIN70_Signature", title = "CIN70 Signature: STK11-MUT vs STK11-WT",
                   color = "STK11_MUT_vs_STK11_WT", 
                   fill = "STK11_MUT_vs_STK11_WT", palette = "jco",
                   alpha = 0.5, font.label = list(size = 50, face = "bold"), ggtheme = theme_bw()) 
#yplot + stat_compare_means()+theme(text = element_text(size = 24))
# yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_fill_manual(breaks = c('ALK_MUT','EGFR_MUT','STK11_MUT','KRAS_MUT','Other_genes_MUT'), 
#                                                                                        values=c('#00CED1','#696969','#3CB371','#E7BB00','#4B0082')) yplot <- ggboxplot(ichorcna_data, x = "STK11_MUT_vs_STK11_WT", y = "CIN70_Signature", title = "CIN70 Signature: STK11-MUT vs STK11-WT",
#                    color = "STK11_MUT_vs_STK11_WT", fill = "STK11_MUT_vs_STK11_WT", palette = "jco",
#                    alpha = 0.5, ggtheme = theme_bw())
#yplot + stat_compare_means()
#colSTKvsNonSTKmain <- c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000')
#yplot + stat_compare_means()+theme(text = element_text(size = 20)) 
yplot + stat_compare_means()+theme(text = element_text(size = 24))+theme_bw()+theme(legend.position = "none")+#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11_MUT_vs_STK11_WT)),
           y = aggregate(CIN70_Signature ~ STK11_MUT_vs_STK11_WT, data_ordered, median)[ , 2],
           label = c('33','11'),
           col = "red",
           vjust = - 1)+ stat_compare_means()+ scale_color_manual(values=c('STK11_WT'='#3CB371','STK11_MUT'='#000000')) +scale_fill_manual(values=c('STK11_WT'='#3CB371','STK11_MUT'='#000000'))
#yplot  + stat_compare_means(method = "t.test")
#ggsave("CIN70_signature_STK11mut_vs_NonSTK11mut.pdf",height = 10,width = 10)
#ggsave("CIN70_signature_STK11mut_vs_NonSTK11mut_1a.pdf",height = 5,width = 5)
ggsave("CIN70_signature_STK11mut_vs_NonSTK11mut_3.pdf",height = 5,width = 5)
# 
# c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')
# ichorcna_data <- read.csv(file="~/Documents/Ben_Izar_Project/fraction_altered_genomes/Fraction_Genome_Altered_STK11_LUAD_TCGA_TMB.csv")
# ichorcna_data <- as.data.frame(ichorcna_data)
# table(ichorcna_data$STK11_mut_vs_Non_STK11_mut)#=='STK11_mut')
# colnames(ichorcna_data)
# yplot <- ggboxplot(ichorcna_data, x = "STK11_mut_vs_Non_STK11_mut", y = "TMB", title = "Tumor mutation burden - TCGA (LUAD)\n STK-mut (All samples) vs Non-STK-mut (KRAS_mut)",
#                    color = "STK11_mut_vs_Non_STK11_mut", fill = "STK11_mut_vs_Non_STK11_mut", palette = "jco",
#                    alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
# yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('Non_STK11_mut'='#3CB371', 'STK11_mut'='#000000')) +scale_fill_manual(values=c('Non_STK11_mut'='#3CB371', 'STK11_mut'='#000000'))

#manual color
ichorcna_data <- read.csv(file="CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
group_ordered <- with(ichorcna_data,                       # Order boxes by median
                      reorder(STK11_MUT_vs_STK11_WT,
                              CIN70_Signature,
                              median))
group_ordered                                     # Print order
data_ordered <- ichorcna_data                              # Create data with reordered group levels
data_ordered$STK11_MUT_vs_STK11_WT <- factor(data_ordered$STK11_MUT_vs_STK11_WT,
                                             levels = levels(group_ordered))
# df<-data_ordered %>%
#   group_by(Barcode) %>%
#   summarise(Freq = sum(STK11_MUT_vs_STK11_WT))


yplot <- ggboxplot(data_ordered, x = "STK11_MUT_vs_STK11_WT", y = "CIN70_Signature", title = "CIN70 Signature: STK11-mut vs Non-STK11-mut",
                   color = "STK11_MUT_vs_STK11_WT", fill = "STK11_MUT_vs_STK11_WT", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#colSTKvsNonSTKmain <- c('STK11_WT'='#3CB371', 'STK11_MUT'='#000000')
yplot + scale_color_manual(values=c('STK11_WT'='#3CB371','STK11_MUT'='#000000')) +scale_fill_manual(values=c('STK11_WT'='#3CB371','STK11_MUT'='#000000'))
yplot + stat_compare_means()+theme(text = element_text(size = 24))+theme_bw()+theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  annotate("text",
           #x = 1:length(table(ichorcna_data$Sample.Type)),
           x = 1:length(table(data_ordered$STK11_MUT_vs_STK11_WT)),
           y = aggregate(CIN70_Signature ~ STK11_MUT_vs_STK11_WT, data_ordered, median)[ , 2],
           label = c('33','11'),
           col = "red",
           vjust = - 1)
#yplot  + stat_compare_means(method = "t.test")
#ggsave("CIN70_signature_STK11mut_vs_NonSTK11mut_v5.pdf",height = 10,width = 10)
ggsave("CIN70_signature_STK11mut_vs_NonSTK11mut_v6.pdf",height = 10,width = 10)

ichorcna_data <- read.csv(file="Integrated_NSCLC_tumor_CX.17.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
for(i in 1:17)
{
  yplot <- ggboxplot(ichorcna_data, x = "STK11mut_vs_NonSTK11mut", y = paste0('CX',i), title = paste0("CX",i,": STK11-mut vs Non-STK11-mut"),
                     color = "STK11mut_vs_NonSTK11mut", fill = "STK11mut_vs_NonSTK11mut", palette = "jco",
                     alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
  yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000')) +scale_fill_manual(values=c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#000000'))
  
  #yplot  + stat_compare_means(method = "t.test")
    ggsave(paste0("CX",i,"_STK11mut_vs_NonSTK11mut_v6.pdf"),height = 10,width = 10)
}
#stress_module Primary_BM_STK_NonSTK
colPrimvsBMvsSTKvsNonSTKmain <- c('Non-STK11-mut (BRAIN_METS)'='#00CED1','Non-STK11-mut (PRIMARY)'='#696969','STK11-mut (BRAIN_METS)'='#FF4500','STK11-mut (PRIMARY)'='#4B0082')


# 6) CIN70 signature Primary_BM_STK_NonSTK
# -----------------------------

#manual color
ichorcna_data <- read.csv(file="Integrated_NSCLC_tumor_CX.17.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
for(i in 1:17)
{
  yplot <- ggboxplot(ichorcna_data, x = "Primary_BM_STK_NonSTK", y = paste0('CX',i), title = paste0("CX",i,": STK11-mut (Primary) vs STK11-mut (BM) \n Non-STK11-mut (Primary) vs Non-STK11-mut (BM) "),
                     color = "Primary_BM_STK_NonSTK", fill = "Primary_BM_STK_NonSTK", palette = "jco",
                     alpha = 0.5, font.label = list(size = 42, face = "bold"), ggtheme = theme_bw())
  yplot + stat_compare_means()+theme(text = element_text(size = 20)) + scale_color_manual(values=colPrimvsBMvsSTKvsNonSTKmain) +scale_fill_manual(values=colPrimvsBMvsSTKvsNonSTKmain)+rotate_x_text(90)
  
  #yplot  + stat_compare_means(method = "t.test")
    ggsave(paste0("CX",i,"_Primary_BM_STK_NonSTK_v6.pdf"),height = 10,width = 10)
  
}
table(ichorcna_data$STK11mut_vs_NonSTK11mut)

# 7) Factor blocks CIN signatures
# -----------------------------

# 7.1) CIN25 signature PRIMARY_vs_BRAIN_METS
#library(ggplot2)
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "primary_vs_bm", y = "CIN25_signature_score_upd", title = "Factor Blocks - CIN25 Signature: PRIMARY vs BRAIN_METS",
                   color = "primary_vs_bm", fill = "primary_vs_bm", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)

# 7.2) CIN25 signature STK11mut_vs_NonSTK11mut
#CIN25 signature STK11mut_vs_NonSTK11mut
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "stkmut_vs_nonstkmut", y = "CIN25_signature_score_upd", title = "Factor Blocks - CIN25 Signature: STK11-mut vs Non-STK11-mut",
                   color = "stkmut_vs_nonstkmut", fill = "stkmut_vs_nonstkmut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN_signature_STK11mut_vs_NonSTK11mut.pdf",height = 10,width = 10)

# 7.3) CIN70 signature PRIMARY_vs_BRAIN_METS
library(ggplot2)
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "primary_vs_bm", y = "CIN70_signature_score_upd", title = "Factor Blocks - CIN70 Signature: PRIMARY vs BRAIN_METS",
                   color = "primary_vs_bm", fill = "primary_vs_bm", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)

# 7.4) CIN70 signature STK11mut_vs_NonSTK11mut
#
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "stkmut_vs_nonstkmut", y = "CIN70_signature_score_upd", title = "Factor Blocks - CIN70 Signature: STK11-mut vs Non-STK11-mut",
                   color = "stkmut_vs_nonstkmut", fill = "stkmut_vs_nonstkmut", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN70_signature_STK11mut_vs_NonSTK11mut.pdf",height = 10,width = 10)


######### 
# 7.5) KINOMO Meta-programs CIN signatures
#CIN25 signature PRIMARY_vs_BRAIN_METS
library(ggplot2)
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
yplot <- ggboxplot(ichorcna_data, x = "FB", y = "CIN25_signature_score_upd", title = "Factor Blocks - CIN25 Signature",
                   color = "FB", fill = "FB", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN25_signature.pdf",height = 10,width = 10)
dev.off()
# #CIN25 signature STK11mut_vs_NonSTK11mut
# ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
# ichorcna_data <- as.data.frame(ichorcna_data)
# colnames(ichorcna_data)
# yplot <- ggboxplot(ichorcna_data, x = "stkmut_vs_nonstkmut", y = "CIN25_signature_score_upd", title = "Factor Blocks - CIN25 Signature: STK11-mut vs Non-STK11-mut",
#                    color = "stkmut_vs_nonstkmut", fill = "stkmut_vs_nonstkmut", palette = "jco",
#                    alpha = 0.5, ggtheme = theme_bw())
# yplot + stat_compare_means()
# #yplot  + stat_compare_means(method = "t.test")
# ggsave("FB_sel_CIN_signature_STK11mut_vs_NonSTK11mut.pdf",height = 10,width = 10)

# 7.6) KINOMO Meta-programs CIN70 signature PRIMARY_vs_BRAIN_METS

library(ggplot2)
ichorcna_data <- read.csv(file="FB_sel_CIN70.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
# my_comparisons <- list( c("FB6", "FB9"), c("FB6", "FB3"), c("FB6", "FB4"), c("FB6", "FB7"), c("FB6", "FB8"), c("FB6", "FB11"),c("FB6", "FB12"),
#                         c("FB9", "FB3"),  c("FB9", "FB4"), c("FB9", "FB7"), c("FB9", "FB8"), c("FB9", "FB11"),c("FB9", "FB12"),
#                         c("FB3", "FB4"), c("FB3", "FB7"), c("FB3", "FB8"), c("FB3", "FB11"),c("FB3", "FB12"),
#                         c("FB4", "FB7"), c("FB4", "FB8"), c("FB4", "FB11"),c("FB4", "FB12"),
#                         c("FB7", "FB8"), c("FB7", "FB11"),c("FB7", "FB12"),
#                         c("FB8", "FB11"),c("FB8", "FB12"),
#                         c("FB11", "FB12"))
yplot <- ggboxplot(ichorcna_data, x = "FB", y = "CIN70_signature_score_upd", title = "Factor Blocks - CIN70 Signature",
                   color = "FB", fill = "FB", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
#yplot + stat_compare_means(comparisons = my_comparisons)
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("FB_sel_CIN70_signature.pdf",height = 10,width = 10)
dev.off()

# 7.7) CIN70 signature Cluster_21 Prim vs BM

ichorcna_data <- read.csv(file="Cluster_21_NSCLC_meta_data.csv")
ichorcna_data <- as.data.frame(ichorcna_data)
colnames(ichorcna_data)
table(ichorcna_data$PRIMARY_vs_BRAIN_METS)
table(ichorcna_data$STK11_mut_vs_Non_STK11_mut)
ichorcna_data[1:10,]
# my_comparisons <- list( c("FB6", "FB9"), c("FB6", "FB3"), c("FB6", "FB4"), c("FB6", "FB7"), c("FB6", "FB8"), c("FB6", "FB11"),c("FB6", "FB12"),
#                         c("FB9", "FB3"),  c("FB9", "FB4"), c("FB9", "FB7"), c("FB9", "FB8"), c("FB9", "FB11"),c("FB9", "FB12"),
#                         c("FB3", "FB4"), c("FB3", "FB7"), c("FB3", "FB8"), c("FB3", "FB11"),c("FB3", "FB12"),
#                         c("FB4", "FB7"), c("FB4", "FB8"), c("FB4", "FB11"),c("FB4", "FB12"),
#                         c("FB7", "FB8"), c("FB7", "FB11"),c("FB7", "FB12"),
#                         c("FB8", "FB11"),c("FB8", "FB12"),
#                         c("FB11", "FB12"))
yplot <- ggboxplot(ichorcna_data, x = "PRIMARY_vs_BRAIN_METS", y = "CIN70_signature_score_upd", title = "Cluster_21 - CIN70 Signature: PRIMARY_vs_BRAIN_METS",
                   color = "PRIMARY_vs_BRAIN_METS", fill = "PRIMARY_vs_BRAIN_METS", palette = "jco",
                   alpha = 0.5, ggtheme = theme_bw())
#yplot + stat_compare_means(comparisons = my_comparisons)
yplot + stat_compare_means()
#yplot  + stat_compare_means(method = "t.test")
ggsave("Cluster_21_CIN70_signature_PRIMARY_vs_BRAIN_METS.pdf",height = 10,width = 10)
dev.off()
