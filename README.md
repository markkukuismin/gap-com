![GitHub](https://img.shields.io/github/license/markkukuismin/gap-com)

# gap-com
General regularization selection method for sparse undirected graphical models with (distinct) cluster structure based on gap statistic.

# Install package (to appear)

For now, use the R-script "gap_com.R".

<!--I used RStudio (version 1.1.453) and Microsoft R Open (version 3.5.1.) to create this package. -->

<!--Unzip the "gapcom.zip" file into a working directory and run the following lines: -->

<!--
library(devtools)
library(roxygen2)
install("gapcom") -->

# Generic algorithm

Gap-com measures the difference between the estimated number of graph clusters and expected number of graph clusters. The expected number of graph clusters can be derived either by 

(i) using so called reference distribution (multivariate uniform, independent variables).

(ii) using so called reference graph (Erdos-Renyi model, aka. random graph model) where the probability for drawing an edge between two arbitrary vertices is defined from the reference distribution.

A generic R implementation of gap-com is straightforward. Assume that the sparse network estimation method corresponds to a function *f* depending on a parameter *lambda* which controls the graph sparsity and *f* returns a list of sparse estimates with increasing sparsity level (the sparsest model first and the densest last), *Y* is an n x p data matrix, *detect.cluster* is the community/cluster detection algorithm, *sparsity* returns the sparsity level of a graph and B is the number of generated reference distributions,

```r
> Graphs = f(Y, lambda)
> nlambda = length(lambda)
> k = rep(0, nlambda)
> for(i in 1:nlambda){
+  Clusters = detect.cluster(Graps[i]) # community detection
+  k[i] = length(Clusters) # nmb of estimated clusters
+}
> Expk = matrix(0, B, nlambda)
> for(b in 1:B){
+ if(useReferenceDist){
+   YNULL = apply(Y, 2, function(x) runif(length(x), min(x), max(x)))
+   RefGraphs = f(YNULL, lambda)
+   for(i in 1:nlambda){
+     RefClusters = detect.cluster(RefGraphs[i]) # community detection
+     Expk[b, i] = length(RefClusters) # nmb of estimated clusters from reference data
+   }
+  }
+ if(useReferenceGraph){
+  if(b == 1){
+   YNULL = apply(Y, 2, function(x) runif(length(x), min(x), max(x)))
+   DummyGraphs = f(YNULL, lambda)
+  }
+   for(i in 1:nlambda){
+       RefGraphs = erdos.renyi.game(p, sparsity(DummyGraphs[i]) , type="gnp") # see igraph
+       RefClusters = detect.cluster(RefGraphs[i]) # community detection
+       Expk[b, i] = length(RefClusters) # nmb of estimated clusters from reference graph
+     }
+ }
+}
> Expk = colMeans(Expk) # The expected nmb of clusters under the (i) reference distribution or (ii) reference graph
> Gap_lambda = Expk - k
> GapIndex = which.max(Gap_lambda)
> opt.lambda = lambda[GapIndex]
```

which returns gap-statistic values (*Gap_lambda*), the index of the largest regularization parameter which maximizes the gap-statistics (*GapIndex*) and the largest value of the regularization parameter which maximizes gap-statistics (*opt.lambda*).

# Example

```r
library(huge)
library(igraph)
library(ggplot2)

source("gap_com.R")

set.seed(6011)

L = huge.generator(d = 200, n = 500, graph = "cluster", g = 7)

Y = L$data

nlambda = 50

HugeSolutionPath = huge(Y, method = "ct", nlambda = nlambda)

gapUnifLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50, method = "unif_sample") # reference distribution (unif sample)
```
![GapComUnif](https://user-images.githubusercontent.com/40263834/83128616-55299e80-a0e4-11ea-80d8-474e38904324.png)

```r
gapERLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50, method = "er_sample") # Erdos-Renyi model
```

![GapComER](https://user-images.githubusercontent.com/40263834/83128641-61156080-a0e4-11ea-800c-1e186f6c0aef.png)

```r
huge.plot(L$theta)

title("Ground truth")
```

![GroundTruthGraph](https://user-images.githubusercontent.com/40263834/83129166-0f210a80-a0e5-11ea-8ecf-44fb45ef64ea.png)

In this example, the graphs selected using the reference distribution resampling or the reference graph resampling are the same because the gap-com statistic is maximized with the same tuning parameter value.

```r

GGapUnif = graph.adjacency(HugeSolutionPath$path[[gapUnifLambda$opt.index]], mode="undirected")

GGapER = graph.adjacency(HugeSolutionPath$path[[gapERLambda$opt.index]], mode="undirected")

graph.isomorphic(GGapUnif, GGapER)
[1] TRUE

huge.plot(HugeSolutionPath$path[[gapLambda$opt.index]])

title("gap-com, unif sample (pairwise correlation hard thresholding)")
```

![GapComUnifGraph](https://user-images.githubusercontent.com/40263834/83129203-1c3df980-a0e5-11ea-8f9a-c6561206c78f.png)

```r
huge.plot(HugeSolutionPath$path[[gapERLambda$opt.index]])

title("gap-com, ER sample (pairwise correlation hard thresholding)")
```

![GapComERGraph](https://user-images.githubusercontent.com/40263834/83129228-27912500-a0e5-11ea-8f83-482fc75971a4.png)

The identified communities are very close to the ground truth clusters:

```r
TrueG = graph.adjacency(L$theta, mode = "undirected", diag = F)

TrueCommunities = walktrap.community(TrueG)

GapERCommunities = walktrap.community(GGapER)

compare(GapERCommunities, TrueCommunities, method="adjusted.rand") # close to one = better
[1] 0.9942448
```

Now also available parallel. Parallel computing is usefull when the number of parameters increas:

```r
L = huge.generator(d = 500, n = 500, graph = "cluster", g = 7)

Y = L$data

HugeSolutionPath = huge(Y, method = "ct", nlambda = nlambda)

# Without parallel:

system.time(GapLambdaER <- gap_com(HugeSolutionPath, B = 50, method = "er_sample"))
   user  system elapsed 
   7.20    0.22    7.42 

# With parallel:

library(foreach)
library(parallel)
library(doParallel)

source("../Rfunctions/gap_com_parallel.R")

registerDoParallel(cores=12)

system.time(GapLambdaERPar <- gap_com_parallel(HugeSolutionPath, B = 50, method = "er_sample"))
   user  system elapsed 
   1.82    0.23    3.63 
```

# Reference

Gap-com statistic is described in:

Kuismin and Sillanpaa (manuscript) "Gap-com: General model selection method for sparse undirected networks with clustering structure".
