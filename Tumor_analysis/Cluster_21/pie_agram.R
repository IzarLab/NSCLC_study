# Author: Somnath Tagore, Ph.D. 
# Title: Pie-agram analysis 
# Script Name: pie_agram.R
# Last Updated: 09/20/2022

### pie-agram
## data input (number of reads mapped to each category)
total=2812
PA048=559
KRAS_6=436
PA042=317
PA001=243
PA068=227
PA104=165
STK_1=82
STK_21=70
PA067=65
KRAS_13=50
PA080=43
KRAS_7=37
STK_14=37
PA019=33
PA034=33
PA125=33
KRAS_8=30
PA070=30
KRAS_4=27
PA056=27
PA005=26
PA054=23
STK_5dot2=22
KRAS_12=19
STK_3=19
KRAS_17=18
N561=18
PA025=17
PA076=17
PA004=14
STK_18=14
PA141=12
KRAS_10=9
PA043=8
PA072=7
KRAS_11=5
N586=4
STK_15=4
STK_20=4
PA060=3
N254=2
STK_5dot1=2
STK_2=1

#contribution
total=100
PA048=28.1332451499118
KRAS_6=7.40390758699373
PA042=1.66585235277543
PA001=27.6337450199203
PA068=1.81972614840989
PA104=1.9973417721519
STK_1=4.04124409448819
STK_21=0.689205673758865
PA067=1.52792048929664
KRAS_13=1.17788640087384
PA080=0.604868654311039
KRAS_7=3.10544041450777
STK_14=0.53824355971897
PA019=0.544486956521739
PA034=1.64235059760956
PA125=0.818073878627968
KRAS_8=0.528035539387155
PA070=0.919420289855073
KRAS_4=1.41561589403973
PA056=0.87379963347587
PA005=0.766817959980478
PA054=0.595583931133429
STK_5dot2=0.54376
KRAS_12=0.580641967422549
STK_3=0.486632243258749
KRAS_17=0.896692307692308
N561=0.600572289156627
PA025=0.694378346222487
PA076=1.38559183673469
PA004=0.487480810864003
STK_18=0.458121256570101
PA141=0.499912070343725
KRAS_10=0.489634990560101
PA043=0.441733671886021
PA072=0.593619402985075
KRAS_11=0.458616874135546
N586=0.503811509591326
STK_15=0.471263616557734
STK_20=0.978080808080808
PA060=0.442371450498849
N254=0.418521212121212
STK_5dot1=0.422157371379161
STK_2=0.41444233807267


cluster_21_samples <- c('PA048',
                        'KRAS_6',
                        'PA042',
                        'PA001',
                        'PA068',
                        'PA104',
                        'STK_1',
                        'STK_21',
                        'PA067',
                        'KRAS_13',
                        'PA080',
                        'KRAS_7',
                        'STK_14',
                        'PA019',
                        'PA034',
                        'PA125',
                        'KRAS_8',
                        'PA070',
                        'KRAS_4',
                        'PA056',
                        'PA005',
                        'PA054',
                        'STK_5dot2',
                        'KRAS_12',
                        'STK_3',
                        'KRAS_17',
                        'N561',
                        'PA025',
                        'PA076',
                        'PA004',
                        'STK_18',
                        'PA141',
                        'KRAS_10',
                        'PA043',
                        'PA072',
                        'KRAS_11',
                        'N586',
                        'STK_15',
                        'STK_20',
                        'PA060',
                        'N254',
                        'STK_5dot1',
                        'STK_2')

cluster_21_samples <- c(PA048,
                        KRAS_6,
                        PA042,
                        PA001,
                        PA068,
                        PA104,
                        STK_1,
                        STK_21,
                        PA067,
                        KRAS_13,
                        PA080,
                        KRAS_7,
                        STK_14,
                        PA019,
                        PA034,
                        PA125,
                        KRAS_8,
                        PA070,
                        KRAS_4,
                        PA056,
                        PA005,
                        PA054,
                        STK_5dot2,
                        KRAS_12,
                        STK_3,
                        KRAS_17,
                        N561,
                        PA025,
                        PA076,
                        PA004,
                        STK_18,
                        PA141,
                        KRAS_10,
                        PA043,
                        PA072,
                        KRAS_11,
                        N586,
                        STK_15,
                        STK_20,
                        PA060,
                        N254,
                        STK_5dot1,
                        STK_2)
PT_BM <- c(BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           PRIMARY,
           BRAIN_METS,
           PRIMARY,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           PRIMARY,
           BRAIN_METS,
           BRAIN_METS,
           CHEST_WALL_MET,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           BRAIN_METS,
           PRIMARY,
           BRAIN_METS,
           BRAIN_METS,
           PRIMARY,
           PRIMARY,
           PRIMARY,
           PRIMARY,
           BRAIN_METS,
           PRIMARY,
           BRAIN_METS,
           PRIMARY)
rRNA=5 # mapped to nuclear rRNA regions
mtRNA=7  # mapped to mitochondria genome
# for the rest of above, then we divide into different category, like http://www.biomedcentral.com/1741-7007/8/149 did.
intergenic=48 
introns=12
exons=30
upstream=3
downstream=6
not_near_genes=40

rest=total-rRNA-mtRNA
genic=rest-intergenic
introns_and_exons=introns+exons-genic

# parameter for pie chart
iniR=0.2 # initial radius
colors=list(NO='grey',total='black',mtRNA='red',rRNA='blue',genic='green',intergenic='yellow',introns='magenta',
            exons='darkorchid',upstream='chocolate4',downstream='bisque1',not_near_genes='darkkhaki')

colPrimvsBMmain <- c('BRAIN_METS'='#FF0000','PRIMARY'='#0000FF')
colPrimvsBMmainTweakBM <- c('BRAIN_METS'='#FFF0F0','PRIMARY'='#0000FF')
colPrimvsBMmainTweakPRIMARY<- c('BRAIN_METS'='#FF0000','PRIMARY'='#EBF7FF')
colSTKvsNonSTKmain <- c('STK11-WT'='#3CB371', 'STK11-MUT'='#000000')
colSTKvsNonSTKTweakSTK <- c('Non-STK11-mut'='#3CB371', 'STK11-mut'='#F1F1F1')
colSTKvsNonSTKTweakNonSTK <- c('Non-STK11-mut'='#EBFFE3', 'STK11-mut'='#000000')

cols = c(
  'Primary'='#A80D11',
  'Brain_Mets'='#008DB8')

cols = c(
  'STK11-MUT'='#A80D11',
  'STK11-WT'='#008DB8')

cols_pt_bm = c(
  'PA048'='#FF0000',
    'KRAS_6'='#FF0000',
    'PA042'='#FF0000',
    'PA001'='#FF0000',
    'PA068'='#FF0000',
    'PA104'='#FF0000',
    'STK_1'='#FF0000',
    'STK_21'='#0000FF',
    'PA067'='#FF0000',
    'KRAS_13'='#0000FF',
    'PA080'='#FF0000',
    'KRAS_7'='#FF0000',
    'STK_14'='#FF0000',
    'PA019'='#FF0000',
    'PA034'='#FF0000',
    'PA125'='#FF0000',
    'KRAS_8'='#FF0000',
    'PA070'='#FF0000',
    'KRAS_4'='#FF0000',
    'PA056'='#FF0000',
    'PA005'='#FF0000',
    'PA054'='#FF0000',
    'STK_5dot2'='#FF0000',
    'KRAS_12'='#0000FF',
    'STK_3'='#FF0000',
    'KRAS_17'='#FF0000',
    'N561'='#ffa500',
    'PA025'='#FF0000',
    'PA076'='#FF0000',
    'PA004'='#FF0000',
    'STK_18'='#FF0000',
    'PA141'='#FF0000',
    'KRAS_10'='#0000FF',
    'PA043'='#FF0000',
    'PA072'='#FF0000',
    'KRAS_11'='#0000FF',
    'N586'='#0000FF',
    'STK_15'='#0000FF',
    'STK_20'='#0000FF',
    'PA060'='#FF0000',
    'N254'='#0000FF',
    'STK_5dot1'='#FF0000',
    'STK_2'='#0000FF'
)
colSTKvsNonSTKmain <- c('STK11-WT'='#3CB371', 'STK11-MUT'='#000000')

cols_stk11_wt = c(
  'PA048'='#3CB371',
    'KRAS_6'='#3CB371',
    'PA042'='#3CB371',
    'PA001'='#3CB371',
    'PA068'='#3CB371',
    'PA104'='#3CB371',
    'STK_1'='#000000',
    'STK_21'='#000000',
    'PA067'='#3CB371',
    'KRAS_13'='#3CB371',
    'PA080'='#3CB371',
    'KRAS_7'='#3CB371',
    'STK_14'='#000000',
    'PA019'='#3CB371',
    'PA034'='#3CB371',
    'PA125'='#3CB371',
    'KRAS_8'='#3CB371',
    'PA070'='#3CB371',
    'KRAS_4'='#3CB371',
    'PA056'='#3CB371',
    'PA005'='#3CB371',
    'PA054'='#3CB371',
    'STK_5dot2'='#000000',
    'KRAS_12'='#3CB371',
    'STK_3'='#000000',
    'KRAS_17'='#3CB371',
    'N561'='#3CB371',
    'PA025'='#3CB371',
    'PA076'='#3CB371',
    'PA004'='#3CB371',
    'STK_18'='#000000',
    'PA141'='#3CB371',
    'KRAS_10'='#3CB371',
    'PA043'='#3CB371',
    'PA072'='#3CB371',
    'KRAS_11'='#3CB371',
    'N586'='#3CB371',
    'STK_15'='#000000',
    'STK_20'='#000000',
    'PA060'='#3CB371',
    'N254'='#3CB371',
    'STK_5dot1'='#000000',
    'STK_2'='#000000'
)
cols = c(
  'PA001' = '#808080', 'PA004' = '#d3d3d3', 'PA005' = '#2f4f4f',
  'PA019' = '#556b2f', 'PA025' = '#8b4513','PA034' = '#2e8b57', 'PA042' = '#228b22',
  'PA043' = '#7f0000','PA048' = '#191970', 'PA054' = '#808000', 'PA056' = '#b8860b',
  'PA060' = '#008b8b',
  'PA067' = '#4682b4', 'PA068' = '#d2691e','PA070' = '#9acd32','PA072' = '#cd5c5c',
  'PA076' = '#00008b', 'PA080' = '#32cd32','PA104' = '#8fbc8f','PA125' = '#8b008b',
  'PA141' = '#b03060','N254' = '#ff4500',#'N561' = '#00ced1',
  'N586' = '#ffa500',
  'KRAS_10' = '#ffd700', 'KRAS_11' = '#6a5acd', 'KRAS_12' = '#deb887',
  'KRAS_13' = '#00ff00', 'KRAS_17' = '#00fa9a', 'KRAS_4' = '#dc143c',
  'KRAS_6' = '#0000ff', 'KRAS_7' = '#a020f0', 'KRAS_8' = '#adff2f',
  'STK_1' = '#da70d6', 'STK_14' = '#ff00ff', 'STK_15' = '#1e90ff',
  'STK_18' = '#f0e68c', 'STK_2' = '#dda0dd','STK_20' = '#90ee90',
  'STK_21' = '#ffa07a', 'STK_22dot2' = '#87cefa', 'STK_3' = '#7fffd4',
  'STK_5dot1' = '#ff69b4','STK_5dot1' = '#ffb6c1'
)
library('plotrix')

# from outer circle to inner circle
#0 circle: blank

pdf(".cluster_21/cluster_21_pie_agram_v3.pdf",height = 10,width=10)
pie(1, radius=iniR, init.angle=90, col=c('beige'), border = NA, labels='')

#4 circle: show genic:exons and intergenic:downstream
# floating.pie(0,0,c(exons, genic-exons+not_near_genes, downstream, mtRNA+rRNA+intergenic-not_near_genes-downstream),
#              radius=5*iniR, startpos=pi/2, col=as.character(colors[c('exons','NO','downstream','NO')]),border=NA)
##floating.pie(0,0,cluster_21_samples, radius=5*iniR, startpos=pi/2, col=as.character(cols_stk11_wt),border=NA)

#3 circle: show genic:introns and intergenic:not_near_genes | upstream
# floating.pie(0,0,c(genic-introns, introns, not_near_genes, intergenic-upstream-not_near_genes, upstream, mtRNA+rRNA),
#              radius=4*iniR, startpos=pi/2, col=as.character(colors[c('NO','introns','not_near_genes','NO','upstream','NO')]),border=NA)
floating.pie(0,0,cluster_21_samples, radius=4*iniR, startpos=pi/2,col=as.character(cols_pt_bm),border=NA)

#2 circle: divide the rest into genic and intergenic
#floating.pie(0,0,c(genic, intergenic, mtRNA+rRNA),radius=3*iniR, startpos=pi/2, col=as.character(colors[c('genic','intergenic','NO')]),border=NA)
floating.pie(0,0,cluster_21_samples,radius=3*iniR, startpos=pi/2, col=as.character(cols),border=NA)

#1 circle: for rRNA+mtRNA+rest
#floating.pie(0,0, c(rest, rRNA,mtRNA), radius=2*iniR, startpos=pi/2, col=as.character(colors[c('NO','rRNA','mtRNA')]), border = NA)

#legend(0, 5*iniR, gsub("_"," ",names(c(cols,cols_pt_bm,cols_stk11_wt))[-1]), col=as.character(cols[-1],cols_pt_bm[-1],cols_stk11_wt[-1]), pch=19,bty='n', ncol=2)
# legend(0,5*iniR,gsub("_"," ",names(c(cols,cols_pt_bm,cols_stk11_wt))[-1]), 
#        col=as.character(cols[-1],cols_pt_bm[-1],cols_stk11_wt[-1]),
#        pch=19,bty='n', ncol=5)
dev.off()
