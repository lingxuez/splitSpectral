#################################################
## Spectral Clustering with Sample Splitting:
## An example on the political blog data.
##
## Copyright Lingxue Zhu (lzhu@cmu.edu)
## All Rights Reserved.
##
## Reference:
## 1. Lei and Zhu (2017),
##   "Generic Sample Splitting for Refined Community Recovery in Degree Corrected Stochastic Block Models",
##   Statistica Sinica
## 2. Adamic and Glance (2005),
##   "The political blogosphere and the 2004 US election: Divided they blog.",
##   Proceedings of the 3rd International Workshop on Link Discovery, ACM, 36–43.
##
###################################################

rm(list=ls())
source("splitSpectral.R")

###################
## load data
###################
## political blog data
PolblogsEdge <- read.csv("data/clean_PolblogsEdges.csv", stringsAsFactors = FALSE)
PolblogsNode <- read.csv("data/clean_PolblogsNodeLabels.csv", stringsAsFactors = FALSE)

## construct the adjacency matrix
n.node <- nrow(PolblogsNode)
BlogAdj <- matrix(0, nrow=n.node, ncol=n.node)
BlogAdj[cbind(PolblogsEdge[,"node1"], PolblogsEdge[,"node2"])] <- 1
BlogAdj[cbind(PolblogsEdge[,"node2"], PolblogsEdge[,"node1"])] <- 1

## 1222 Nodes: 586 liberal and 636 conservative
table(PolblogsNode[, "truelabel"])
## 1   2 
## 586 636 

## 16714 Edges
sum(BlogAdj[upper.tri(BlogAdj)])
## 16714


######################
## Community detection
######################

## spectral clustering
spectral.clust.est = SpectralClust(Adj=BlogAdj, K=2, isSphere=TRUE)
spectral.accur = Accuracy(spectral.clust.est, PolblogsNode[, "truelabel"], K=2)
print(paste("Sphere spectral clustering: accuracy =", spectral.accur))

## 2-fold cross clustering
cross.2v.clust.est = CrossClust.vFold(Adj=BlogAdj, fold=2, K=2, isSphere=TRUE) 
cross.2v.accur = Accuracy(cross.2v.clust.est, PolblogsNode[, "truelabel"], K=2)
print(paste("2-fold sphere cross clustering: accuracy =", cross.2v.accur))

## 10-fold cross clustering
cross.10v.clust.est = CrossClust.vFold(Adj=BlogAdj, fold=10, K=2, isSphere=TRUE) 
cross.10v.accur = Accuracy(cross.10v.clust.est, PolblogsNode[, "truelabel"], K=2)
print(paste("10-fold sphere cross clustering: accuracy =", cross.10v.accur))

## self-cross clustering based on spectral clustering on the full graph
cross.self.clust.est = CrossClust(Adj=BlogAdj, 
                                  i.G1=c(1:n.node), G1.clust.est=spectral.clust.est, 
                                  i.G2=c(1:n.node), K=2, isSphere=TRUE)
cross.self.accur = Accuracy(cross.self.clust.est, PolblogsNode[, "truelabel"], K=2)
print(paste("Sphere self-cross clustering: accuracy =", cross.self.accur))

