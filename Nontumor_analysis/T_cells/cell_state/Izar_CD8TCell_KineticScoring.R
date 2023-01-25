library(readr)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
library(cowplot)
library(emmeans)
library(nlme)
library(dplyr)
library(plotrix)
library(rstatix)
library(ggrepel)

load("KineticSets_v4_2YF.rda")

## Backscoring Up Minus Down #########################################################################################################

SignatureScore <- function(InputGeneMatrix,TestSignatures){
  GeneMatrixFiltered = merge(data.frame(TestSignatures),data.frame(InputGeneMatrix),by.x = 1,by.y = 0)
  
  GeneMatrixFiltered_Zscore = data.matrix(GeneMatrixFiltered[,2:ncol(GeneMatrixFiltered)])
  for(i in 1:nrow(GeneMatrixFiltered_Zscore)) GeneMatrixFiltered_Zscore[i,] = scale(GeneMatrixFiltered_Zscore[i,])
  
  MeanZ_Score = numeric(ncol(GeneMatrixFiltered_Zscore))
  
  ## Averages Z score of all genes in a signature (Using Geom mean will return an error for negative Z scores!!!)
  for(i in 1:ncol(GeneMatrixFiltered_Zscore)) MeanZ_Score[i] = mean(GeneMatrixFiltered_Zscore[,i])
  names(MeanZ_Score) = colnames(GeneMatrixFiltered_Zscore)
  
  MeanZ_Score
}

SignatureScore_UpMinusDown <- function(InputGeneMatrix,UpDownList){
  UpSignature = SignatureScore(InputGeneMatrix,UpDownList[["Up"]])
  if(!is.null(UpDownList[["Down"]])){
    DnSignature = SignatureScore(InputGeneMatrix,UpDownList[["Down"]])
    NetSignature = UpSignature - DnSignature
  } else{
    NetSignature = UpSignature
  }
  NetSignature
}

SignatureScore_RowInd <- function(InputGeneMatrix,RowInds){
  GeneMatrixFiltered = InputGeneMatrix[RowInds,]
  
  for(i in 1:nrow(GeneMatrixFiltered)) GeneMatrixFiltered[i,] = scale(GeneMatrixFiltered[i,])
  
  MeanZ_Score = colMeans(GeneMatrixFiltered)
  names(MeanZ_Score) = colnames(InputGeneMatrix)
  
  MeanZ_Score
}

SignatureScore_Random <- function(InputGeneMatrix,UpDownList){
  UpGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Up"]]))
  DnGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Down"]]))
  
  UpSignature = SignatureScore_RowInd(InputGeneMatrix,UpGenesRandom)
  if(!is.null(UpDownList[["Down"]])){
    DnSignature = SignatureScore_RowInd(InputGeneMatrix,DnGenesRandom)
    NetRandom = UpSignature - DnSignature
  } else{
    NetRandom = UpSignature
  }
  NetRandom
}

# KineticSets = KineticLists;SigList = DEGList_UpDown;pal = palette;nperm = 100;SameScale = T;PValsAnnot = T; PValsSize = 5
SignatureScorePlot <- function(SigList,KineticSets = KineticLists,pal = rep("black",length(SigList)),nperm = 100,
                               SameScale = F, PValsAnnot = T,PValsSize = 4.2){
  
  ## Calculates Signature Score and Random Score Matrix (nperm controls number of permutations)
  ScoringResults = lapply(SigList, function(Sig) {
    SigScore_Kinetic = lapply(KineticSets, function(KineticSet){
      SigScore = SignatureScore_UpMinusDown(KineticSet$Exprs, Sig)
      
      RandScore = matrix(nrow = nperm, ncol = ncol(KineticSet$Exprs))
      for(Perms in 1:nperm) RandScore[Perms,] = SignatureScore_Random(KineticSet$Exprs,Sig)
      colnames(RandScore) = colnames(KineticSet$Exprs)
      
      list(SigScore = SigScore, RandScore = RandScore, palette = NULL, PlotType = KineticSet$PlotType)
    })})
  
  ## Adds Color Palette to ScoringResults Object
  for(i in 1:length(SigList)){
    for(j in 1:length(ScoringResults[[i]])){
      ScoringResults[[i]][[j]]$palette = pal[i]
    }
  }
  
  ScoringResults = unlist(ScoringResults,recursive = F)
  
  #### Plotting Functions for different Kinetic Set Types

  ## Plotting Function
  ScorePlots = lapply(ScoringResults,function(x){
    ## If Experiment has replicates
    if(x$PlotType=="Replicate_Line"){
      Timepoints = unique(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
      
      ## Signature Score Calculation
      LME_Model = data.frame(Sample = names(x$SigScore),Exp = x$SigScore,ExpBase = "Exp")
      LME_Model$Rep = sapply(strsplit(names(x$SigScore),split = "_"),function(y) y[length(y)])
      LME_Model$Sample = sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
      LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
      
      LME_Sig_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
      LME_Sig_Results = summary(emmeans(LME_Sig_Results,"pairwise"~Sample,adjust = "none"))
      PlottingMatrix_Sig = LME_Sig_Results[["emmeans"]]
      PlottingMatrix_Sig$lower.CL = PlottingMatrix_Sig$emmean - (PlottingMatrix_Sig$SE*1.96)
      PlottingMatrix_Sig$upper.CL = PlottingMatrix_Sig$emmean + (PlottingMatrix_Sig$SE*1.96)
      PlottingMatrix_Sig = PlottingMatrix_Sig %>% select(c("Sample","emmean","SE","lower.CL","upper.CL"))
      
      ######### P Value Calculation
      PValues = LME_Sig_Results$contrasts
      PValues$Sample1 = factor(sapply(strsplit(PValues$contrast,split = " - "),function(x) x[1]),levels = Timepoints)
      PValues$Sample2 = factor(sapply(strsplit(PValues$contrast,split = " - "),function(x) x[2]),levels = Timepoints)
      PValues = subset(PValues, as.numeric(Sample2) - as.numeric(Sample1) == 1)
      PValues$PVal_Asterisk = as.character(symnum(PValues$p.value,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","")))
      
      ## Random Score Calculation
      x$RandScore = x$RandScore[!is.nan(x$RandScore[,1]),]
      RandScore = list()
      for(i in 1:nrow(x$RandScore)) RandScore[[i]] = x$RandScore[i,]
      LME_Rand_Results = sapply(RandScore,function(x) {
        LME_Model = data.frame(Sample = names(x),Exp = x,ExpBase = "Rand")
        LME_Model$Rep = sapply(strsplit(names(x),split = "_"),function(y) y[length(y)])
        LME_Model$Sample = sapply(strsplit(names(x),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
        LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
        
        LME_Rand_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
        LME_Rand_Results = summary(emmeans(LME_Rand_Results,"pairwise"~Sample,adjust = "none"))[["emmeans"]]$emmean
        names(LME_Rand_Results) = Timepoints
        LME_Rand_Results
      })
      
      LME_Rand_Results = reshape2::melt(LME_Rand_Results,value.name = "emmeans", varnames = c("Sample","Rep"))
      LME_Rand_Results$Sample = factor(LME_Rand_Results$Sample,levels = Timepoints)
      LME_Rand_Results = LME_Rand_Results %>% group_by(Sample) %>%
        summarize(emmean = mean(emmeans), SE = std.error(emmeans), lower.CL = emmean - 1.96*SE, upper.CL = emmean + 1.96*SE)
      
      
      ######### Plotting
      Plot = ggplot(data = LME_Rand_Results,aes(x = Sample,y = emmean,group = 1)) +
        geom_ribbon(aes(ymin = lower.CL,ymax = upper.CL),fill = "grey") +
        geom_line(color = "black") + 
        geom_line(data = PlottingMatrix_Sig,color = x$palette) +
        geom_errorbar(data = PlottingMatrix_Sig,aes(ymin = lower.CL,ymax = upper.CL),color = x$palette,width=0.4, size=0.5) +
        geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      Results = list(Plot = Plot,PVals = PValues, AxisLimits = layer_scales(Plot)$y$get_limits(), PlotType = x$PlotType, palette = x$palette)
      
      ## Else if there are no replicate samples (EG Kupper)  
    }else if(x$PlotType=="Slope_Line"){
      
      ######### P Values
      ## Signature Slope
      Signature_TRMCourse = x$SigScore[-(1:2)]
      Signature_TRMCourse = data.frame(ZScore = Signature_TRMCourse,
                                       Sample = names(Signature_TRMCourse))
      Signature_TRMCourse$Sample = as.numeric(c("0","5","10","15","20","25","30","45","60","90"))
      
      
      LmFit = lm(ZScore ~ Sample,data = Signature_TRMCourse)
      
      SlopeResults = list(DeltaExpr = summary(LmFit)$coefficients[2,1]*90,
                          PVals = summary(LmFit)$coefficients[2,4])
      SlopeResults$PVals = as.character(symnum(SlopeResults$PVals,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","ns")))
      
      
      
      SlopeResults = paste0("DExprs_d0-d90: ",signif(SlopeResults$DeltaExpr,3),", Sig: ",SlopeResults$PVals)
      
      ######### Plotting
      Timepoints = factor(names(x$SigScore),levels = names(x$SigScore))
      PlottingMatrix_Sig = data.frame(Timepoints,Means = x$SigScore,
                                      Timepoints_v2 = c(-20,-10,0,5,10,15,20,25,30,45,60,90))
      
      x$RandScore = x$RandScore[!is.nan(x$RandScore[,1]),]
      Rand_Results = melt(t(x$RandScore),value.name = "ZScore", varnames = c("Sample","Rep"))
      Rand_Results$Sample = factor(Rand_Results$Sample, levels = Timepoints)
      Rand_Results = Rand_Results %>% group_by(Sample) %>%
        summarize(Means = mean(ZScore), SE = std.error(ZScore), LowerCL = Means - 1.96*SE, UpperCL = Means + 1.96*SE)
      
      Plot = ggplot(data = Rand_Results,aes(x = Timepoints,y = Means,group = 1)) + 
        geom_ribbon(fill = "grey",mapping = aes(ymin = LowerCL, ymax = UpperCL)) +
        geom_line(color = "black") + 
        geom_smooth(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],method = "lm", formula = y~x, color = "transparent") +
        geom_line(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],color = x$palette) +
        geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
        geom_smooth(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],method = "lm",formula = y~x,color = colorspace::darken(x$palette,0.2),fill = "transparent") +
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      Results = list(Plot = Plot,PVals = SlopeResults, AxisLimits = layer_scales(Plot)$y$get_limits(), PlotType = x$PlotType, palette = x$palette)
    }else if(x$PlotType=="Replicate_Line_TTest"){
      Timepoints = unique(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
      
      ##########################
      ## Stats
      
      SigScore = melt(x$SigScore, value.name = "ZScore")
      SigScore$Sample = factor(sapply(strsplit(rownames(SigScore),split = "_"),function(x) x[1]),levels = Timepoints)
      SigScore_PVals = SigScore %>% pairwise_t_test(ZScore ~ Sample, p.adjust.method = "none",pool.sd = T) %>% 
        subset(as.numeric(factor(group2,levels = Timepoints)) - as.numeric(factor(group1,levels = Timepoints)) == 1)
      SigScore_PVals$p.signif = as.character(symnum(SigScore_PVals$p,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","")))
      
      
      ##########################
      ## Random Calculation
      
      x$RandScore = x$RandScore[!is.nan(x$RandScore[,1]),]
      Rand_Results = melt(t(x$RandScore),value.name = "ZScore", varnames = c("Sample","Rep"))
      Rand_Results$Sample = factor(sapply(strsplit(rownames(SigScore),split = "_"),function(x) x[1]),levels = Timepoints)
      Rand_Results = Rand_Results %>% group_by(Rep, Sample) %>%
        summarize(ZScore_Bootstrap = mean(ZScore), n = length(ZScore),.groups = "keep")
      Rand_Results = Rand_Results %>% group_by(Sample) %>%
        summarize(Means = mean(ZScore_Bootstrap), SE = std.error(ZScore_Bootstrap), UpperCL = Means + 1.96*SE, LowerCL = Means - 1.96*SE)
      
      ##########################
      ## Signature Calculation
      
      SigScore = SigScore %>% group_by(Sample) %>%
        summarize(Means = mean(ZScore), SE = std.error(ZScore), UpperCL = Means + 1.96*SE, LowerCL = Means - 1.96*SE)
      
      ##########################
      ## Plotting
      
      Plot = ggplot(data = Rand_Results,aes(x = Sample,y = Means,group = 1)) + 
        geom_ribbon(data = Rand_Results,aes(ymin = LowerCL,ymax = UpperCL),fill = "grey") +
        geom_line(data = Rand_Results,color = "black") + 
        geom_line(data = SigScore, color = x$palette) +
        geom_point(data = SigScore, fill = x$palette,shape = 21,size = 4) +
        geom_errorbar(data = SigScore, aes(ymin = LowerCL,ymax = UpperCL),color = x$palette,width=0.4, size=0.5) + 
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      
      Results = list(Plot = Plot,PVals = SigScore_PVals, AxisLimits = layer_scales(Plot)$y$get_limits(), PlotType = x$PlotType, palette = x$palette)
    }else{
      Results = NULL
      }})
  
  ## Adjusts/Syncs Scales if User Specified (SameScale Variable)
  if(SameScale){
    AxisMin = min(sapply(ScorePlots,function(x) x$AxisLimits[1]))
    AxisMax = max(sapply(ScorePlots,function(x) x$AxisLimits[2]))
    
    ScorePlots = lapply(ScorePlots, function(x){
      if(!PValsAnnot) x$Plot = x$Plot + 
          ylim(c(AxisMin,AxisMax))
      
      x$AxisLimits = c(AxisMin,AxisMax)
      x
    })
    
    NoAxisLabels = 1:length(ScorePlots)
    for(i in NoAxisLabels[-(((length(SigList)-1)*length(KineticSets))+1)]) ScorePlots[[i]]$Plot = ScorePlots[[i]]$Plot + theme(axis.text.y = element_blank())
  }
  
  ## Adds P Value Annotation if User Specified
  if(PValsAnnot){
    ScorePlots = lapply(ScorePlots, function(x){
      
      AxisLimit_PVal = (x$AxisLimits[2] - x$AxisLimits[1])*1.15 + x$AxisLimits[1]
      
      if(x$PlotType=="Replicate_Line"){
        x$Plot = x$Plot + 
          geom_text(data = x$PVals, mapping = aes(x = Sample2, y = x$AxisLimits[2], label = PVal_Asterisk),size = PValsSize,color = colorspace::darken(x$palette,0.2),vjust = "bottom") +
          ylim(c(x$AxisLimits[1],AxisLimit_PVal))
      }else if(x$PlotType=="Slope_Line"){
        x$Plot = x$Plot + 
          annotate("text",label = x$PVals,x = 1, y = mean(c(x$AxisLimits[2],AxisLimit_PVal)),hjust = "left",color = colorspace::darken(x$palette,0.2),vjust = "bottom") + 
          ylim(c(x$AxisLimits[1],AxisLimit_PVal))
      }else if(x$PlotType=="Replicate_Line_TTest"){
        x$Plot = x$Plot + 
          geom_text(data = x$PVals,mapping = aes(x = group2,y = x$AxisLimits[2],label = p.signif),size = PValsSize,color = colorspace::darken(x$palette,0.2),vjust = "bottom") +
          ylim(c(x$AxisLimits[1],AxisLimit_PVal))
      }
      return(x)
    })
  }
  
  ## Combines Plots
  
  CombinedPlot = lapply(ScorePlots, function(x) x$Plot)
  
  ## Makes Pretty = removes x axis from top plots (just on bottom x axis), adjusts x axis labels
  PlotIndices = 1:length(CombinedPlot)
  for(i in 1:(7*(length(SigList)-1))) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
  for(i in (7*(length(SigList)-1)+1):length(CombinedPlot)) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.text.x = element_text(size = 6))
  

  plot_grid(plotlist= CombinedPlot, nrow = length(SigList),align = "hv",rel_widths = c(0.7,1,1,1.2,1.4,0.5,0.7))
}

## Imports Query Sets and Kinetic Sets  #########################################################################################################

DEGList = msigdbi::read.gmt("Izar_CD8TCell_DEGs_v3.gmt")$genesets
DEGList = DEGList[c(1,2,5,6)]

DEGList_UpDown = list()

for(i in seq(1,4,by = 2)) DEGList_UpDown[[names(DEGList)[i]]] = list(Up = DEGList[[i]],Down = DEGList[[i+1]])
names(DEGList_UpDown) = sapply(strsplit(names(DEGList_UpDown),split = ".",fixed = TRUE),function(x) x[1])

palette = c("black","blue")

Scale = 1.9
pdf("Izar_CD8TCell_DEGs_KineticScoring.pdf",width = 6.55*Scale,height = 2*Scale)
SignatureScorePlot(DEGList_UpDown,pal = palette, SameScale = T)
dev.off()








