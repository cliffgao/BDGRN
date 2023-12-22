# Bayesian Directed GRN for Single Cell Data
# (R implementation.)

# Authors:
# Yao Mengting Nankai University

#! /usr/bin/env Rscript

source('/BDGRN/code/KNN.R')
source('/BDGRN/code/BaCo.R')

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(modEvA))
suppressPackageStartupMessages(library("psych"))


BDGRN = function(data, t, cutoff =0, k = 3, w = 0.1){
## data is SingleCellData which is a matrix with m genes by n cells 
## t is timeseries of cells in data
## cutoff control the threshold of Bayesian Correlations to get edgelists
## k is parameter of KNN, choosing k nearest neighbors
## w is parameter of sliding windows, ranging 0-1 .

#pseudo-time/real-time
t[,"time"] = apply(t,1,max,na.rm=TRUE)
data = data[,order(t$time)] #对数据进行时间排序
name<-rownames(data)
a = knn_smoothing(data,k=k) #knn_smooth
#rownames(a) = name
## 基因表达矩阵a,窗口长度s，中间矩阵m*（ncol(a)-s+1),相关系数矩阵a_corr
genes=nrow(a)
cells=ncol(a)
s =ceiling(cells*w)#滑窗长度

cunchu=matrix(0,(cells-s+1),genes)
a_corr=matrix(0,genes,genes)
for (e in 1:nrow(a)){

gene_e =a[e,1:s]
for (k in 0:(cells-s)){

gene_k=a[1:genes,(1+k):(k+s)]
gene_k[e,1:s]=gene_e  #相应位置替换为固定窗口
g_cor=BaCo(gene_k)
cunchu[k+1,1:genes]=g_cor[e,1:genes]
}

cor_max=apply(cunchu,2,max)
a_corr[e,1:genes]=cor_max

}

dimnames(a_corr)=list(name,name)

#相关矩阵
b<-abs(a_corr)
b[b==1]<-0
g1<-graph.adjacency(abs(b),mode="directed",weighted=T)
e1<-get.edgelist(g1)
e2<-cbind(e1,round(E(g1)$weight,3))

threshold = quantile(b[b!=0],cutoff)
e2[which(e2[1:nrow(e2),3]<threshold),3]=0
colnames(e2) =c("Gene1","Gene2","Probability")

edgelist = data.frame(e2)
edgelist

}
