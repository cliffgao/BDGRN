#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(modEvA))
suppressPackageStartupMessages(library("psych"))
suppressPackageStartupMessages(library("PRROC"))


#计算aupr
comAUPR = function( e_test, e_pred){
## e_test are edges in ground-truth network
## e_pred are edges predicted by BDGRN

e_test<-e_test[,1:2]
e_test=as.matrix(e_test)
e_test<-e_test[which(e_test[,1]!=e_test[,2]),]
e_test<-data.frame(e_test)
m<-paste(e_test$Gene1,e_test$Gene2,sep="") #设置标准网络边的标签g1g2.

e_pred$Gene1<-paste(e_pred$Gene1,e_pred$Gene2,sep="")
e_pred$Gene2 = 0
midx1<-match(m,e_pred$Gene1)
e_pred[na.exclude(midx1),2]<-1

trueEdges = e_pred[1:nrow(e_pred),2]
predEdges= e_pred[1:nrow(e_pred),3]
pr_out=pr.curve(predEdges, weights.class0=trueEdges)
roc=roc.curve(predEdges, weights.class0=trueEdges)
aupr=pr_out$auc.integral
auroc=roc$auc

res = list(aupr,auroc)
res

}


#计算early precision
comEP = function(e_test, e_pred){
## e_test are edges in ground-truth network
## e_pred are edges predicted by BDGRN

e_test<-e_test[,1:2]
e_test=as.matrix(test)
e_test<-e_test[which(e_test[,1]!=e_test[,2]),]
e_test<-data.frame(e_test)
e_pred =e_pred[order(e_pred[,3],decreasing =T),]
e_pred = e_pred[1:nrow(e_test),1:2]

c<-intersect(e_test,e_pred)
ep=nrow(c)/nrow(e_test)

ep
}

    
