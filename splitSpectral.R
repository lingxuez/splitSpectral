#################################################
## Spectral Clustering with Sample Splitting
##
## Copyright Lingxue Zhu (lzhu@cmu.edu)
## All Rights Reserved.
## 
## Reference:
## Lei and Zhu (2017),
## "Generic Sample Splitting for Refined Community Recovery in Degree Corrected Stochastic Block Models",
## Statistica Sinica
##
#################################################

library(mgcv)

RowNorm <- function(A) {
  ## Normalize the non-zero rows of a matrix A
  ## to have unit L2 norm
  ##
  ## Parameters:
  ## --A: the input matrix
  ##
  ## Returns:
  ## the new matrix after row-normalization.
  
  A.rownorm <- sqrt(rowSums(A^2))
  i.nonzero <- which(A.rownorm > 0)
  A[i.nonzero, ] <- apply(A[i.nonzero, ], 2, function(x){x / A.rownorm[i.nonzero]})
  return (A)
}


SpectralClust <- function(Adj, K=2, isSphere=FALSE){
  ## Spectral clustering
  ##
  ## Parameters:
  ## --Adj: adjacency matrix
  ## --K: number of communities
  ## --isSphere: whether to use sphere spectral clustering for degree corrected models
  ##
  ## Returns:
  ## a vector of estimated cluster memberships.
  
  U <- mgcv::slanczos(Adj, K)$vectors
  if (isSphere) {
    ## normalize non-zero rows to account for node heterogeneity
    U <- RowNorm(U)
  }
  clust.est <- kmeans(U, centers=K, nstart=20)$cluster
  return(clust.est)
}

CrossClust <- function(Adj, i.G1, G1.clust.est, i.G2, K=2, isSphere=FALSE) {
  ## Cross clustering:
  ## use the estimated membership on i.G1 to cluster i.G2
  ##
  ## Parameters:
  ## --Adj: adjacency matrix
  ## --i.G1: indices of nodes in the first group
  ## --G1.clust.est: cluster membership of i.G1
  ## --i.G2: indices of nodes in the second group, disjoint from i.G1
  ## --K: number of communities
  ## --isSphere: whether to use sphere methods for degree corrected models
  ##
  ## Returns:
  ## a vector of estimated memberships on i.G2
  
  ## the connectivity between each node in i.G2 and the K clusters in i.G1
  U.G2 <- matrix(nrow=length(i.G2), ncol=K)
  for (i in 1:nrow(U.G2)) {
    for (k in 1:K) {
      U.G2[i, k] = mean(Adj[i.G2[i], i.G1[G1.clust.est==k]])
    }
  }
  if (isSphere) {
    ## normalize non-zero rows to account for node heterogeneity
    U.G2 = RowNorm(U.G2)
  }
  ## clustering
  G2.clust.est <- kmeans(U.G2, centers=K, nstart=20)$cluster
  return(G2.clust.est)
}


CrossClust.vFold <- function(Adj, 
                             fold=2, K=2, 
                             isSphere=FALSE ## whether to perform sphere clustering
                             ) {
  ## V-fold cross clustering
  ##
  ## Parameters:
  ## --Adj: adjacency matrix
  ## --fold: number of folds to use, at least 2
  ## --K: number of communities
  ## --isSphere: whether to use sphere methods for degree corrected models
  ##
  ## Returns:
  ## a vector of estimated memberships
  
  ## randomly split nodes into V folds
  n <- nrow(Adj)
  fold.size <- ceiling(n/fold)
  permute.index <- sample(1:n, size=n, replace=FALSE)
  fold.index <- split(permute.index, ceiling((1:n)/fold.size))
  
  ## estimate membership on every fold using CrossClust
  fold.clust.est <- vector("list", fold)
  for (i in 1:fold){
    i.G2 <- fold.index[[i]]
    i.G1 <- c(1:n)[-i.G2] 
    ## spectral clustering on i.G1
    G1.clust.est <- SpectralClust(Adj=Adj[i.G1, i.G1], K=K, isSphere=isSphere)
    ## cross clustering on i.G2
    G2.clust.est <- CrossClust(Adj=Adj, i.G1=i.G1, G1.clust.est=G1.clust.est,
                                     i.G2=i.G2, K=K, isSphere=isSphere)
    fold.clust.est[[i]] <- G2.clust.est
  } 
  
  ## merge all folds
  merged.clust.est = Merge.vFold(Adj, K, fold, fold.index, fold.clust.est, isSphere)
  
  ## final results: need to re-permute the memberships to get the original order
  final.clust.est <- rep(0, n)
  final.clust.est[unlist(fold.index)] <- unlist(merged.clust.est)
  return(final.clust.est)

}

Merge.vFold <- function(Adj, K, fold, 
                        fold.index, fold.clust.est, isSphere=FALSE) {
  ## Merge the memberships in all folds
  ##
  ## Parameters:
  ## --Adj: adjacency matrix
  ## --fold: number of folds to use, at least 2
  ## --K: number of communities
  ## --fold.index: a list with length=fold, 
  ##         the i-th element is the indices of nodes in fold-i
  ## --fold.clust.est: a list with length=fold, 
  ##         the i-th element is the estimated memberships for nodes in fold-i
  ## --isSphere: whether to use sphere methods for degree corrected models
  ##
  ## Returns:
  ## a new list with length=fold, where the i-th element
  ## contains the new memberships for nodes in fold-i
  
  merged.clust.est = vector("list", fold)
    
  ## fold-1: reference
  i.fold.ref <- fold.index[[1]]
  clust.est.ref <- fold.clust.est[[1]]
  merged.clust.est[[1]] <- fold.clust.est[[1]]
  
  ## fold-2 to fold-V: permute labels to be consistent with fold-1
  for (i in 2:fold){
    i.fold <- fold.index[[i]]
    clust.est <- fold.clust.est[[i]]
    ## calculate Bhat in Algorithm 3
    B.hat <- EstimateB(Adj=Adj, i.fold.ref=i.fold.ref, i.fold=i.fold, 
                  clust.est.ref=clust.est.ref, clust.est=clust.est, 
                  K=K, isSphere=isSphere)
    ## get the permuted labels
    greedy.perm <- MergeGreedyPerm(B.hat)
    merged.clust.est[[i]] <- Permute(clust.est, greedy.perm, K)
  }
  
  return(merged.clust.est)
}


EstimateB <- function(Adj, i.fold.ref, i.fold, clust.est.ref, clust.est, K, isSphere) {
  ## Estimate the matrix B as in Algorithm 3 for merging labels.
  ##
  ## Parameters:
  ## --Adj: adjacency matrix
  ## --i.fold.ref: the indices of nodes in fold-1 (the reference fold)
  ## --clust.est.ref: the memberships of i.fold.ref
  ## --i.fold: the indices in the current fold that will be merged
  ## --clust.est: the memberships of i.fold
  ## --K: number of communities
  ## --isSphere: whether to use sphere methods for degree corrected models
  ## 
  ## Returns:
  ## the estimated K-by-K matrix B, indicating the connectivities among clusters
  
  B.hat <- matrix(nrow=K, ncol=K) 
  for (i in 1:K) {
    for (j in 1:K) {
      B.hat[i, j] <- mean(Adj[i.fold[clust.est == i], i.fold.ref[clust.est.ref == j]])
    }
  }
  if (isSphere) {
    ## normalize non-zero rows to account for node heterogeneity
    B.hat <- RowNorm(B.hat)
  }
  return (B.hat)
}


MergeGreedyPerm <- function(B.hat) {
  ## Find the best permutation of rows of B.hat
  ## Here, a greedy method is implemented
  ## which is specifically for the case when diag(B) dominates others,
  ## i.e., within-community connectivities are larger than between-community connectivities
  ## This function permutes the rows of B.hat such that diag(B.hat) is the largest
  ##
  ## Parameters:
  ## --B.hat: the estimated matrix B
  ##
  ## Returns:
  ## the permutation indices
  
  K <- nrow(B.hat)
  greedy.perm <- rep(0, K)
  
  ## permute labels for clust.est
  while (sum(B.hat!=0) > 0) {
    ## put the next largest value on the diagonal
    i.max <- which.max(B.hat)
    max.row <- row(B.hat)[i.max]
    max.col <- col(B.hat)[i.max]
    greedy.perm[max.row] <- max.col
    
    ## remove the merged labels
    B.hat[max.row, ] <- 0
    B.hat[, max.col] <- 0
  }
  
  ## for the remaining un-handled communities, randomly permute
  if (sum(greedy.perm==0) > 0) {
    clust.remaining <- c(1:K)[-greedy.perm[greedy.perm!=0]]
    greedy.perm[greedy.perm==0] <- clust.remaining
  }
  return (greedy.perm)
}


Permute <- function(clust.est, greedy.perm, K) {
  ## permute clust.est such that
  ## the original i-th cluster becomes greedy.perm[i]
  ##
  ## Parameters:
  ## --clust.est: the original membership vector
  ## --greedy.perm: the permutation
  ## --K: number of communities
  ##
  ## Returns:
  ## the permuted membership vector, where cluster i becomes cluster greed.perm[i]
  
  new.clust.est <- rep(0, length(clust.est))
  for (k in 1:K) {
    new.clust.est[which(clust.est == k)] <- greedy.perm[k]
  }
  return(new.clust.est)
}


Accuracy <- function(clust.est, clust.true, K){
  ## Calculating accuracy, 
  ## using a greedy method to find the optimal permutation of labels
  ##
  ## Parameters:
  ## --clust.est: the estimated cluster membership
  ## --clust.true: the true cluster membership
  ## --K: number of clusters
  ##
  ## Returns the accuracy

  inter.len <- table(clust.est, clust.true)
  
  ## permute labels for clust.est
  greedy.perm <- rep(0, K)
  while (sum(inter.len!=0) > 0) {
    ## find the next largest intersection
    i.max <- which.max(inter.len)
    max.row <- row(inter.len)[i.max]
    max.col <- col(inter.len)[i.max]
    greedy.perm[max.row] <- max.col
    
    ## remove from consideration
    inter.len[max.row, ] <- 0
    inter.len[, max.col] <- 0
  }
  
  if (sum(greedy.perm==0) > 0) {
    ## remaining ones: randomly assign
    clust.remaining <- c(1:K)[-greedy.perm[greedy.perm!=0]]
    greedy.perm[greedy.perm==0] <- clust.remaining
  }
  
  ## accuracy
  accur <- sum(clust.true == Permute(clust.est, greedy.perm, K)) / length(clust.true)
  return(accur)
}



