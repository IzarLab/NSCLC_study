# Author: Somnath Tagore, Ph.D. Title: Random Forest Classifier
# Script Name: classifier.R
# Last Updated: 03/19/2022

library(Hmisc)
library("caret")
library("rpart")
library("tree")
library("e1071")
library(ggplot2) # Data visualization
library(readr) # CSV file I/O, e.g. the read_csv function
library(randomForest)
library(randomForestExplainer)
library(caret)
library(randomForest)
library(mlbench)
library(MLeval)

# Input data files are available in the "../input/" directory.
# For example, running this (by clicking run or pressing Shift+Enter) will list the files in the input directory

#system("ls ../input")
pat<-'STK_21'
testdata<-read.csv(paste0("./random_forest/nsclc.kras.stk.classifier/testing_data/",pat,"/",pat, "_testing_rf.csv"),header=TRUE, sep=",")
rownames(testdata) <- testdata[,1]
testdata<-testdata[,-1]

Sonar <- cl.2.target
inTraining <- createDataPartition(Sonar$Class, p = .8, list = FALSE)
#inTraining <- createDataPartition(Sonar[,101], p = .5, list = FALSE)
training <- Sonar[ inTraining,]
testing  <- Sonar[-inTraining,]

#testing <- mix.pik3ca.brca.cl.1.2.3.4.wt
testing <1- testdata
fitControl <- trainControl(## 5-fold CV
  method = "cv",
  summaryFunction=twoClassSummary, 
  classProbs=T,
  savePredictions = T,
  number = 5)

#set.seed(825)
rfFit1 <- train(Class ~ ., data = training, 
                 method = "rf", 
                 trControl = fitControl,
                 metric = "ROC",
                 verbose = TRUE)
rfFit1
plot(rfFit1)
fm <- rfFit1$finalModel


importance <- varImp(rfFit1)

#importance(fm)
# estimate variable importance
#importance <- varImp(rfFit1, scale=FALSE)
importance
#write.csv(importance$importance,file=paste0("~/Documents/random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "e_importance_rf.csv"))
#write.csv(importance$importance,file=paste0("~/Documents/random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_ichorCNA_evidence_importance_rf.csv"))
write.csv(importance$importance,file=paste0("./random_forest/nsclc.kras.stk.classifier/testing_data/",pat,"/",pat, "_ichorCNA_evidence_importance_rf.csv"))
# summarize importance
print(importance)
# plot importance

pdf(paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_importance_rf.pdf"))
plot(importance,20)
dev.off()



predict.rf=predict(rfFit1,testing,type="prob")

write.csv(predict.rf,file=paste0("./random_forest/nsclc.kras.stk.classifier/testing_data/",pat,"/",pat, "_ichorCNA_evidence_predict_rf.csv"))
res <- evalm(rfFit1)
ggsave(
  paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_ROC_rf.pdf"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 7,
  height = 7,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

importance <- varImp(svmFit1, scale=FALSE)
write.csv(importance$importance,file=paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_importance_rf.csv"))
# summarize importance
print(importance,20)
# plot importance
pdf(paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_importance_rf.pdf"))
plot(importance,20)
dev.off()
# 
predict.svm=predict(svmFit1,testing,type="prob")
write.csv(predict.svm,file=paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_predict_rf.csv"))
res <- evalm(svmFit1)
ggsave(
  paste0("./random_forest/nsclc.classifier/testing_data/",pat,"/",pat, "_ROC_rf.pdf"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
