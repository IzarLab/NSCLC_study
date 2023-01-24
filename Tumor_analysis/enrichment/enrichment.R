# Author: Somnath Tagore, Ph.D. 
# Title: Enrichment Analysis (GSEA, Pathway)
# Script Name: enrichment.R
# Last Updated: 09/01/2022

# Load packages
#install.packages("enrichR")
library(enrichR)
library(ggplot2)
library(dplyr)

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# combine both Cancer hallmarks (ch) and GO

# gene set
mylist <- c('ANKRD36C',
            'C3',
            'CXCL2',
            'LPCAT1',
            #MT-ATP6,
            'NDRG1',
            'SFTPA2',
            'SFTPC',
            'SYNE1',
            'BCCIP',
            'BTBD9',
            'IGFBP7',
            'FN1',
            'MALAT1',
            'MGP',
            'SPARC',
            'C7',
            'CCN1',
            'COL6A3',
            'HBB',
            'IGFBP2',
            'IGKC',
            'ROBO2',
            'TIMP3',
            'CAMK2N1',
            'CD82',
            'CEACAM5',
            'CLU',
            'DANT2',
            'LMO7',
           # MT-CO2
            #MT-ND2
            #MT-ND3
            'NAV2',
            'NDUFC1',
            'NUCB2',
            'PLAAT4',
            'REV3L',
            'SCEL',
            'SCGB3A1',
            'SCGB3A2',
            'SFTPA1',
            'SLPI',
            'TFPI',
            'VPS53',
            'AHNAK',
            'AHNAK2',
            'ANXA1',
            'ANXA5',
            'AQP3',
            'AUTS2',
           # C12orf57
            'CAMK1D',
            'CAV2',
            'CELSR1',
            'CHMP2A',
            'COX14',
           # CP
            'CXCL14',
            'CXCL17',
            'DPYSL2',
            'DRAM1',
            'EEA1',
            'ESF1',
            'HMGB2',
            'HMGCS1',
            'HMGN5',
            'HS6ST2',
            'HSPE1',
            'IER2',
            'IFT57',
            'JUN',
            'KRT19',
            'LAMB3',
            'LDLR',
            'MAGOH',
            'MAP2',
            'MCTP2',
           # MT-CO1
            #MT-ND1
          #  MT-ND4
           # MT-ND5
          #  MT2A
            'MYADM',
            'NDUFB4',
            'PARP8',
            'PDE3B',
            'PSMA7',
            'PTMS',
            'PTPN13',
            'PTPN4',
            'PVT1',
            'QSOX1',
            'RPL14',
            'SCD',
            'SFTPB',
            'SKIL',
            'TC2N',
            'TNFRSF10B',
            'TNFSF10',
            'TRIO',
            'UGCG',
            'YBX1',
            'ZNF710',
            'ABCC4',
            'ABL1',
            'ADCY9',
            'ADIRF',
            'AGAP1',
            'AGPAT3',
            'AIG1',
            'ALOX15B',
            'ANGPTL4')
tumor.mylist <- mylist

# list which pathway sets need to scored 
dbs <- c("MSigDB_Hallmark_2020","GO_Biological_Process_2021")
if (websiteLive) {
  #enriched <- enrichr(unique(mylist[1:100,1]), dbs)
  #enriched <- enrichr(mylist[,1], dbs)
  enriched <- enrichr(mylist, dbs)
}

#mut_enr<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))

# calculate qscore
mut_enr_ch<-mutate(enriched[[1]], qscore = -log(Adjusted.P.value, base=10))
mut_enr_go<-mutate(enriched[[2]], qscore = -log(Adjusted.P.value, base=10))

# create a data frame with all pathways scored by qscore
dim(mut_enr_ch)
dim(mut_enr_go)
mut_enr <- rbind.data.frame(mut_enr_ch,mut_enr_go)
dim(mut_enr)
celltype = "Tumor"

enrichdbs = "Cancer_Hallmarks_and_Gene_Ontology"
# topn = "Top100"
# mp = "MP15"
topn = "Top100"
#mp = "MP15"
h_mut_enr1 <- mut_enr[1:25,]#[1:500,] 
h_mut_enr <- h_mut_enr1
#df[df$col1 == 1, ]

# plot the result
ggp<- h_mut_enr %>%
  #extract(Term, into = c("Term", "id"), "(.*)\\s\\((.*)\\)") %>%
  #mutate(Term=stri_trans_totitle(Term)) %>%
  ggplot(aes(qscore, reorder(Term, qscore), fill = P.value)) +
  scale_fill_gradient(low = "red", high = "blue") +
  geom_bar(stat = "identity") +
  
  #facet_wrap( ~ Class, nrow = 3, , scales = "free") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank()
  ) +
  geom_text(
    aes(label = paste("P.val=", round(P.value,3))),
    color = "black",
    size = 4,
    hjust = 1, nudge_x = 2
    #+ 
    # xlab("qscore") + ylab("Description")
    #position = position_dodge(width = 1.5)
  ) + theme_bw()
myggp<-ggp + scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+ggtitle(paste(mp,": ",enrichdbs,"\ntop = ",topn,"_genes")) +xlab("qscore") + ylab("Description")
pdf(paste0(celltype,"_Enrichment_",enrichdbs,"_",topn,"_v1.pdf"), width = 10, height = 10)
print(myggp)
#print(mygg1p)
#print(mygg2p)
dev.off()
