#! /usr/bin/env Rscript
#PCA降维

#定义Bayesian相关性函数
BaCo <- function(X){ 
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}

#按列归一化
normalize=function(x){return((x-min(x))/(max(x)-min(x)))}

#PCA降维函数
#prcomp()和princomp()
#PCA(){FactoMineR包}
#dudi.pca()[ade4包]
#epPCA()[ExPosition包]



library(igraph)
library(dplyr)
library(modEvA)
library("psych")

win_len= seq(0.1,0.9,0.1)
setwd("/public/home/yaomt/BEELINE_data/")
path_list=dir()
  
auroc_matrix=matrix(nrow=10,ncol=length(win_len))
rownames(auroc_matrix)=path_list
for (l in win_len){
for (j in 1:length(path_list)){

load_path=paste("/public/home/yaomt/BEELINE_data/",path_list[j],sep="")
print(load_path)
setwd(load_path)
path<-dir()
path=path[path!="Rplots.pdf"]
auroc_list<-matrix(0,nrow=length(path),ncol=1)
for (i in 1:length(path)){

#导入标准网络

test<-read.csv(file=paste(path[i],"/refNetwork.csv",sep=""),header=T)
test<-test[,1:2]
g_ref<-graph_from_edgelist(as.matrix(test),directed=F)
e_test<-get.edgelist(g_ref)
e_test<-e_test[which(e_test[,1]!=e_test[,2]),]
#e_test<-unique(e_test)
e_test<-data.frame(e_test)


m<-paste(e_test$X1,e_test$X2,sep="")
n<-paste(e_test$X2,e_test$X1,sep="")

#e_test

a<-read.csv(file=paste(path[i],"/ExpressionData.csv",sep=""),header=F,row.names=1)
a<-a[-1,]
name<-rownames(a)
a<-apply(a,2,as.numeric)
rownames(a)<-name
print(i)
## 基因表达矩阵a,窗口长度s，中间矩阵m*（ncol(a)-s+1),相关系数矩阵a_corr
genes=nrow(a)
cells=ncol(a)
s =ceiling(cells*l)
print(paste("l:",l))
print(paste("s:",s))
cunchu=matrix(0,(cells-s+1),genes)
print(dim(cunchu))
a_corr=matrix(0,genes,genes)
dimnames(a_corr)=list(name,name)
for (e in 1:nrow(a)){
gene_e =a[e,1:s]
#b=a[-e,:]
for (k in 0:(cells-s)){
gene_k=a[1:genes,(1+k):(k+s)]
gene_k[e,1:s]=gene_e  #相应位置替换为固定窗口
g_cor=BaCo(gene_k)
cunchu[k+1,1:genes]=g_cor[e,1:genes]
}
cor_max=apply(cunchu,2,max)
#print(cor_max)
a_corr[e,1:genes]=cor_max
}
#dim(a_corr)
dimnames(a_corr)=list(name,name)

#相关矩阵
#a_corr<-BaCo(a)
b<-abs(a_corr)
b[b==1]<-0
cutoff=median(b[b!=0])
print(paste("cutoff:",cutoff))
b[b<cutoff]=0
g1<-graph.adjacency(abs(b),mode="directed",weighted=T)
e1<-get.edgelist(g1)
e2<-cbind(e1,round(E(g1)$weight,3))
e2<-data.frame(e2)
e2$X1<-paste(e2$X1,e2$X2,sep="")
e2$X2<-0
midx1<-match(m,e2$X1)
midx2<-match(n,e2$X1)
e2[na.exclude(midx1),2]<-1
e2[na.exclude(midx2),2]<-1


#计算auroc
k<-e2[,2:3]
k<-as.matrix(k)
k<-apply(k,2,as.double)
k<-data.frame(k)
#print(k)
auroc<-modEvA::AUC(obs=k$X2,pred=k$X3,curve="ROC",simplif=TRUE,main="ROC")
print(auroc)
auroc_list[i,1]=auroc
}

median_auroc<-apply(auroc_list,2,median)
print(paste("median_auroc",median_auroc))
print(j)
auroc_matrix[j,10*l]=median_auroc
}
}
print(auroc_matrix)
#write.csv(aupr_ratio_matrix,file="/public/home/yaomt/pca_aupr_result4.csv",row.names=F)
