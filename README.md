# splitSpectral
This is the R code for 
>  Lei and Zhu (2017), "Generic Sample Splitting for Refined Community Recovery in Degree Corrected Stochastic Block Models", Statistica Sinica.

Our algorithm uses a sample splitting trick with spectral clustering for refined community recovery in stochastic block models.
Pease cite our work in your publication if it helps your research:
```
@article{Lei2017generic,
    title={Generic Sample Splitting for Refined Community Recovery in Degree Corrected Stochastic Block Models},
    author={Lei, Jing and Zhu, Lingxue},
    journal={Statistica Sinica},
    year={2017}
}
```

To re-produce our results on the [political blog data](https://networkdata.ics.uci.edu/data.php?id=102), 
please refer to the `R` script
```
example_polblogs.R
```

To use our code in `R`:
```{r}
source("splitSpectral.R")

## here is a toy example with 100 nodes, forming 2 fully connected groups, each with 50 nodes
## we construct the adjacency matrix
test_adj = matrix(0, nrow=100, ncol=100)
test_adj[1:50, 1:50] = 1
test_adj[51:100, 51:100] = 1
true_clusters = c(rep(1, 50), rep(2, 50))

## now we perform cross clustering with sample splitting
K = 2 ## 2 communities
fold = 2 ## 2-fold sample splitting
clusters = CrossClust.vFold(Adj=test_adj, fold=fold, K=K,
    isSphere=FALSE ## set to TRUE for degree-corrected block models
    )

## check the accuracy
print(paste0("Accuracy = ", Accuracy(clusters, true_clusters, K=2)))
```
