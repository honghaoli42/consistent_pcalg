# PC algorithm with consistent separating sets


This repository contains the code for the modification to constraint-based algorithms for causal discovery with consistent separating sets proposed in Li et al. [0], applied to the PC-stable [1] implementation of the pcalg R package [2].

`consistent_utils.R` contains the SimpleGraph class with the necessary components for fast recovery of consistent separating nodes, as detailed in the supplementary information of [0].

`consistent_pcalg_orientation.R` and `consistent_pcalg.R` respectively correspond to algorithm 3 and 4 of [0] for enforcing orientation and skeleton-level consistency. They both overwrite the `skeleton` and `pc` functions of pcalg [2] to return fully consistent graphs.

Dependencies : [`igraph`](https://cran.r-project.org/web/packages/igraph/index.html), [`dequer`](https://cran.r-project.org/web/packages/dequer/index.html)

## Examples

For this example, we use the Hepar2 network from the [bayesian network repository](http://www.bnlearn.com/bnrepository/discrete-large.html#hepar2) and compare the outputs of the original PC-stable and its skeleton-consistent and orientation-consistent versions.
The p-value thresholds are chosen as the ones that maximize the skeleton F1-score, cf Figure 5 of [0].

``` {R}
# Generate test data and convert for pcalg use
library(bnlearn)
Hepar2 <- read.bif(file='hepar2.bif.gz', debug = FALSE)
data <- rbn(Hepar2, n=1000, debug = F)
for(col in colnames(data)){
    data[[col]] = as.numeric(data[[col]])-1
}

# Load pcalg and test data
library(pcalg)
V <- colnames(data) # labels aka node names
## define sufficient statistics
suffStat <- list(dm = data,
                 nlev = apply(data, 2, function(col){length(unique(col))}),
                 adaptDF = FALSE)

#####
## estimate CPDAG with original PC-stable
pc.fit <- pc(suffStat,
             ## independence test: G^2 statistic
             indepTest = disCItest, alpha = 0.5, labels = V, verbose = F)

# Load SimpleGraph class for recovering consistent separating sets
source('path/to/consistent_utils.R')

#####
## Skeleton consistent
source('path/to/consistent_pcalg.R')
skel_consistent_pc.fit <- pc(suffStat,
             ## independence test: G^2 statistic
             indepTest = disCItest, alpha = 5.7e-3, labels = V, verbose = F)

#####
## Orientation consistent
source('path/to/consistent_pcalg_orientation.R')
orientation_consistent_pc.fit <- pc(suffStat,
             ## independence test: G^2 statistic
             indepTest = disCItest, alpha = 5.7e-3, labels = V, verbose = F)
```


## References

[0] [Li H, Cabeli V, Sella N, Isambert H. Constraint-based Causal Structure Learning with Consistent Separating Sets. In Advances in Neural Information Processing Systems 2019 (pp. 14257-14266).](http://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets) | [Poster](https://drive.google.com/open?id=1g1vzhgbAsRjGanS0zLWZzmBdrCdoBY_6)

[1] [Colombo D, Maathuis MH. Order-independent constraint-based causal structure learning. The Journal of Machine Learning Research. 2014 Jan 1;15(1):3741-82.](http://www.jmlr.org/papers/volume15/colombo14a/colombo14a.pdf)

[2] [Kalisch M, Mächler M, Colombo D, Maathuis MH, Bühlmann P. Causal inference using graphical models with the R package pcalg. Journal of Statistical Software. 2012 May 17;47(11):1-26.](ftp://ftp.sam.math.ethz.ch/sfs/pub/Manuscripts/buhlmann/pcalg-jss.pdf)