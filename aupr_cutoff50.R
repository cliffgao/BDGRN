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

density<-read.csv("/public/home/yaomt/density.csv",header=F,row.names=1)
density<-density[-1,]
name<-rownames(density)
density<-apply(density,2,as.numeric)
rownames(density)<-name

setwd("/public/home/yaomt/BEELINE_data/")
path_list=dir()
  
aupr_ratio_matrix=matrix(nrow=10,ncol=1)
rownames(aupr_ratio_matrix)=path_list

for (j in 1:length(path_list)){


load_path=paste("/public/home/yaomt/BEELINE_data/",path_list[j],sep="")
print(load_path)

setwd(load_path)
path<-dir()
aupr_list<-matrix(0,nrow=length(path),ncol=1)
for (i in 1:length(path)){
#print(path[i])
#}

#导入标准网络

test<-read.csv(file=paste(path[i],"/refNetwork.csv",sep=""),header=T)
test<-test[,1:2]
g_ref<-graph_from_edgelist(as.matrix(test),directed=F)
e_test<-get.edgelist(g_ref)
e_test<-e_test[which(e_test[,1]!=e_test[,2]),]
e_test<-unique(e_test)
e_test<-data.frame(e_test)

#计算density
nodes_count<-vcount(g_ref)
edge_count<-nrow(e_test)
#density<-(2*edge_count)/(nodes_count*(nodes_count-1))
net_density=density[j,2]


m<-paste(e_test$X1,e_test$X2,sep="")
n<-paste(e_test$X2,e_test$X1,sep="")

#e_test

a<-read.csv(file=paste(path[i],"/ExpressionData.csv",sep=""),header=F,row.names=1)
a<-a[-1,]
name<-rownames(a)
a<-apply(a,2,as.numeric)
rownames(a)<-name
print(i)
a_corr=BaCo(a)
#相关矩阵
#a_corr<-BaCo(a)
b<-abs(a_corr)
b[b==1]<-0
cutoff=median(b[b!=0])
print(cutoff)
b[b<cutoff]=0
g1<-graph.adjacency(abs(b),mode="undirected",weighted=T)
e1<-get.edgelist(g1)
e2<-cbind(e1,round(E(g1)$weight,3))
e2<-data.frame(e2)
e2$X1<-paste(e2$X1,e2$X2,sep="")
e2$X2<-0
midx1<-match(m,e2$X1)
midx2<-match(n,e2$X1)
e2[na.exclude(midx1),2]<-1
e2[na.exclude(midx2),2]<-1


#计算aupr
k<-e2[,2:3]
k<-as.matrix(k)
k<-apply(k,2,as.double)
k<-data.frame(k)
aupr<-modEvA::AUC(obs=k$X2,pred=k$X3,curve="PR",simplif=TRUE,main="PR")
#print(aupr)
aupr_list[i,1]=aupr
}

print(paste(nodes_count,edge_count))
#print(aupr_list)
median_aupr<-apply(aupr_list,2,median)
print(paste("median_aupr",median_aupr))
print(paste(path_list[j],net_density))
aupr_ratio<-median_aupr/net_density
print(j)
aupr_ratio_matrix[j,]=aupr_ratio
print(paste("auprc_ratio:",aupr_ratio))
#print(median(aupr_list))
}
print(aupr_ratio_matrix)
#write.csv(aupr_ratio_matrix,file="/public/home/yaomt/pca_aupr_result4.csv",row.names=F)
