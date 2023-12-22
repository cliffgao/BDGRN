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


#EPR的pca聚类
library(igraph)
library(dplyr)
library(modEvA)
library("FactoMineR")
#library("factoextra")
library("psych")

density<-read.csv("/public/home/yaomt/density.csv",header=F,row.names=1)
density=density[-1,]
name=rownames(density)
density=apply(density,2,as.numeric)
rownames(density)=name
net_list=c("GSD","HSC","mCAD","VSC")

EPR=matrix(nrow=4,ncol=1)
rownames(EPR)=net_list
for(j in 1:length(net_list)){
net_density=density[net_list[j],2]
count<-matrix(0,nrow=10,ncol=1)

setwd(paste("/public/home/yaomt/BEELINE_data/",net_list[j],sep=""))
path_list=dir()
for (i in 1:length(path_list)){

test<-read.csv(file=paste(path_list[i],"/refNetwork.csv",sep=""),header=T)
test<-test[,1:2]
g_ref<-graph_from_edgelist(as.matrix(test),directed=F)
e_test<-get.edgelist(g_ref)
e_test<-e_test[which(e_test[,1]!=e_test[,2]),]
e_test<-unique(e_test)



a<-read.csv(file=paste(path_list[i],"/ExpressionData.csv",sep=""),header=F,row.names=1)
a<-a[-1,]
name<-rownames(a)
a<-apply(a,2,as.numeric)
rownames(a)<-name
a_corr=BaCo(a)
b<-abs(a_corr)
b[b==1]<-0
cutoff=median(b[b!=0])
print(cutoff)
b[b<cutoff]=0
b[order(abs(b),decreasing=T)[2*nrow(e_test)+1:length(b)]]<-0
b[b>0]<-1
b[b<0]<--1
g1<-graph.adjacency(abs(b),mode="undirected")
e1<-get.edgelist(g1)
e1<-data.frame(e1)
e1_t<-e1 %>% select(X2,X1)
c1<-intersect(data.frame(e_test),e1)
colnames(e1_t)<-c("X1","X2")
c2<-intersect(data.frame(e_test),e1_t)
c<-union(c1,c2)
dim(c)
ep=nrow(c)/nrow(e_test)
count[i,1]<-ep
}
print(count)
median_ep<-apply(count,2,median)
EPR[j,]=median_ep/net_density
print(paste(net_list[j],EPR[j,]))
}
print(EPR)
#write.csv(EPR,file="/public/home/yaomt/EPR_result4.csv",row.names=T)


