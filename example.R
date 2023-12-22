source("/BDGRN/code/BDGRN.R")
source("/BDGRN/code/compute.R")
# 读取单细胞表达数据gene*cell

data = read.csv("/BDGRN/data/example/ExpressionData.csv",header =T ,row.names =1)
##读取细胞时间序列

time = read.csv("/BDGRN/data/example/PseudoTime.csv",header=T,row.names =1)
 ## 计算得到网络的edgelist及probability
##可自行设置参数 k, cutoff, w
edgelist = BDGRN(data,time,k=3,cutoff = 0, w =0.1)
edgelist

##compute AUPR/AUROC/Early Precision
##导入ground-truth network
test<-read.csv("/BDGRN/data/example/refNetwork.csv",header=T)

res = comAUPR(test, edgelist)
EP = comEP(test, edgelist)
