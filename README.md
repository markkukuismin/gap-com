![GitHub](https://img.shields.io/github/license/markkukuismin/gap-com)

# gap-com
General regularization selection method for sparse undirected graphical models with (distinct) cluster structure based on gap statistic.

# Install package

I used RStudio (version 1.1.453) and Microsoft R Open (version 3.5.1.) to create this package.

Unzip the "gapcom.zip" file into a working directory and run the following lines:

```r
library(devtools)
library(roxygen2)

install("gapcom")
```

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
![Gap_comPlot](https://user-images.githubusercontent.com/40263834/74665426-5f183900-51a8-11ea-9f86-38192f65ae3f.png)

```r
gapERLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50, method = "er_sample") # Erdos-Renyi model
```


```r
huge.plot(L$theta)

title("Ground truth")
```
![GroundTruthGraph](https://user-images.githubusercontent.com/40263834/74665693-d948bd80-51a8-11ea-8cb3-6a5ad7409402.png)

```r
huge.plot(HugeSolutionPath$path[[gapLambda$opt.index]])

title("gap-com, unif sample (pairwise correlation hard thresholding)")
```
![EstimatedGraph](https://user-images.githubusercontent.com/40263834/74665753-f1b8d800-51a8-11ea-87ff-630c0874413a.png)

# Reference

Gap-com statistic is described in:

Kuismin and Sillanpaa (manuscript) "Gap-com: General model selection method for sparse undirected networks with clustering structure".
