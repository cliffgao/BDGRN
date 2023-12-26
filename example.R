source("./code/BDGRN.R")
source("./code/compute.R")


# reading in scRNA-seq dataset gene*cell
data = read.csv("./data/example/ExpressionData.csv",header =T ,row.names =1)

# reading the PseudoTime 
time = read.csv("./data/example/PseudoTime.csv",header=T,row.names =1)

#run BDGRN
edgelist = BDGRN(data,time,k=3,cutoff = 0, w =0.1)
print(edgelist)

