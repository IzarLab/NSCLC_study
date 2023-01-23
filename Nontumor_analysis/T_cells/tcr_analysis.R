#!/usr/bin/env Rscript

### title: TCR Analysis for NSCLC
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(stringr)
library(reshape2)
library(patchwork)
library(RColorBrewer)
'%notin%' <- Negate('%in%')

system('aws s3 cp s3://nscl-seq/Seurat/integrated_data/44_samples/nontumor/manual_annotation/group1/NSCLC_44_samples_T_cells_v7.rds nsclc/')

#setwd('~/Documents/nsclc/')
label<-'NSCLC'
celltype<-'tcells'
filename<-paste0(label,'_',celltype)
path.ct <- paste0('nsclc/',celltype,'/')
directory<-'nsclc/tcells/'
ifelse(!dir.exists(file.path(directory)), dir.create(file.path(directory),recursive = T), FALSE)

colExp<-c('#C82F71','#7473A6')
colSha<-c('#a67473','#2fc886','#2f71c8')
colMAIT<-c('#009F3F','#E0C126','#F08080','#2fbec8')
colExpPub<- c('#d85a90','#a8285f','#807fae','#504f7c')


seu<-readRDS(paste0('nsclc/NSCLC_44_samples_T_cells_v7.rds'))
#seu<-subset(seu, cell_type_fine != 'NK cells')
seu$barcode<-rownames(seu@meta.data)
seu$TCR<-ifelse(is.na(seu$CTaa)==F,'TCR','no_TCR')
seu$both_chains<-ifelse(grepl("_NA|NA_",seu$CTaa),'both',
                        ifelse(is.na(seu$CTaa)==F,'not_both',NA))
quantile(seu$Frequency,na.rm=T,probs=c(0.05,0.5,0.95))
q5<-quantile(seu$Frequency,na.rm=T,probs=c(0.05))
seu$clone_size<-ifelse(seu$Frequency > q5, 'expanded',
                       ifelse(is.na(seu$CTaa)==F,'not_expanded',NA))



#### Plot TCR data ####

df <- seu@meta.data[, c('barcode', 'patient', 'patient_long','Frequency', 'clone_size','cloneType', 'CTaa', 'manual.annot.fine', 
                        'both_chains', 'Primary_vs_BM','STK11mut_vs_NonSTK11mut')]
df_both <- df %>% filter(both_chains == 'both') %>% filter(is.na(cloneType)==F)


pdf(file = paste0(directory,'plots_',filename,'_TCR.pdf'))
DimPlot(seu,group.by="manual.annot.fine")
DimPlot(seu,group.by="manual.annot")
print(FeaturePlot(seu, features = c("rna_CD4",'rna_CD8A','rna_TOX','rna_TCF7'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = F)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
print(FeaturePlot(seu, features = c('rna_NCAM1','rna_NKG7','rna_CXCR5','rna_BTLA','rna_FOXP3','rna_SELL'), 
                  min.cutoff = "q05",max.cutoff = "q95", order = T,raster = F)& 
        scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral'))))
DimPlot(seu,group.by="cloneType")
FeaturePlot(seu, features = c('Frequency'), 
            order = T,raster = F)& 
  scale_colour_gradientn(colours = rev(brewer.pal(11,'Spectral')))
hist(seu$Frequency,breaks=300)
DimPlot(seu,group.by="TCR")
DimPlot(seu,group.by="both_chains")
#DimPlot(seu,group.by="clone_size")

# Stacked bar plot both_chains
ggplot(df, aes(x = patient_long, fill = both_chains)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle('TCR chains recovered across patients') + 
  ylab('Fraction')

# Stacked bar plot cloneType
df_clone_size <- df_both %>% group_by(patient_long,cloneType) %>% summarise(n = n())
ggplot(df_clone_size, aes(x = patient_long,y=n, fill = cloneType)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# T cells') 

# Stacked bar plot cloneType CD8
df_clone_size_cd8 <- df_both %>% filter(grepl('CD8',manual.annot.fine)) %>% group_by(patient_long,cloneType) %>% summarise(n = n())
ggplot(df_clone_size_cd8, aes(x = patient_long,y=n, fill = cloneType)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# CD8+ T cells') 

# Bar plot clones
df_unique <- df_both %>% group_by(patient_long) %>% summarise(nClones = n_distinct(CTaa))
p1<-ggplot(df_unique, aes(x = patient_long, y = nClones)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Clones')

# Bar plot cells
df_cells <- df_both %>% group_by(patient_long) %>% summarise(nCells = n_distinct(barcode))
p2<- ggplot(df_cells, aes(x = patient_long, y = nCells)) + 
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Cells')

# Stacked bar plot cloneType
p3<-ggplot(df_both, aes(x = patient_long, fill = cloneType)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ylab('Fraction of T cells')

p1/p2/p3

# Circle stacked bar plot stratified clonotype + Primary_vs_BM
df_ct_strat <- df_both %>%
  group_by(Primary_vs_BM, cloneType, manual.annot.fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = manual.annot.fine, 
                        fill = manual.annot.fine, color = manual.annot.fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(cloneType ~ Primary_vs_BM)+
  theme(legend.position = 'right',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified clonotype + STK11mut_vs_NonSTK11mut
df_ct_strat <- df_both %>%
  group_by(STK11mut_vs_NonSTK11mut, cloneType, manual.annot.fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = manual.annot.fine, 
                        fill = manual.annot.fine, color = manual.annot.fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(cloneType ~ STK11mut_vs_NonSTK11mut)+
  theme(legend.position = 'right',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified cell type + Primary_vs_BM
df_clonesize_strat <- df_both %>%
  group_by(Primary_vs_BM, manual.annot.fine, cloneType) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = cloneType, fill = cloneType, color = cloneType)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(manual.annot.fine ~ Primary_vs_BM)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified cell type + STK11mut_vs_NonSTK11mut
df_clonesize_strat <- df_both %>%
  group_by(STK11mut_vs_NonSTK11mut, manual.annot.fine, cloneType) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = cloneType, fill = cloneType, color = cloneType)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() +
  facet_grid(manual.annot.fine ~ STK11mut_vs_NonSTK11mut)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))

# Circle stacked bar plot stratified Primary_vs_BM + STK11mut_vs_NonSTK11mut
df_clonesize_strat <- df_both %>%
  group_by(STK11mut_vs_NonSTK11mut, Primary_vs_BM, cloneType) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = cloneType, fill = cloneType, color = cloneType)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(Primary_vs_BM ~ STK11mut_vs_NonSTK11mut)+
  theme(legend.position = 'bottom',
        strip.text.x = element_text(angle = 90))
dev.off()
