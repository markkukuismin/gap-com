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

(i) using so called reference distribution (permuted data set).

(ii) using so called reference graph (Erdos-Renyi model, aka. random graph model) where the probability for drawing an edge between two arbitrary vertices is defined from the permuted data.

A generic R implementation of gap-com is straightforward. Assume that the sparse network estimation method corresponds to a function *f* depending on a parameter *lambda* which controls the graph sparsity and *f* returns a list of sparse estimates with increasing sparsity level (the sparsest model first and the densest last), *Y* is an n x p data matrix, *detect.cluster* is the community/cluster detection algorithm, *sparsity* returns the sparsity level of a graph and B is the number of permuted data sets,

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
+   YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
+   RefGraphs = f(YNULL, lambda)
+   for(i in 1:nlambda){
+     RefClusters = detect.cluster(RefGraphs[i]) # community detection
+     Expk[b, i] = length(RefClusters) # nmb of estimated clusters from reference data
+   }
+  }
+ if(useReferenceGraph){
+  if(b == 1){
+   YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
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

source("../RFunctions/gap_com.R")

##########################################################

# Initialize seed number:

#seed = Sys.time()

#seed = as.integer(seed)

#seed = seed %% 100000

seed = 29734

set.seed(seed)

##########################################################

L = huge.generator(d = 200, n = 500, graph = "cluster", g = 7)

Y = L$data

nlambda = 50

HugeSolutionPath = huge(Y, method = "ct", nlambda = nlambda)

gapPermuteLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50, method = "permute_sample") # reference distribution (permuted data set)
```
![GapComPermute](https://user-images.githubusercontent.com/40263834/108598453-5dbeb500-7396-11eb-8369-d8d122c1cfcd.png)

```r
gapERLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50, method = "er_sample") # Erdos-Renyi model
```

![GapComER](https://user-images.githubusercontent.com/40263834/108598460-6dd69480-7396-11eb-8146-80717391fe96.png)

```r
huge.plot(L$theta)

title("Ground truth")
```

![GroundTruthGraph](https://user-images.githubusercontent.com/40263834/108598471-7929c000-7396-11eb-81c1-ee75ad6b24ce.png)

In this example, the graphs selected using the reference distribution resampling or the reference graph resampling are not exactly the same but very close to each other. Actually, the community structre of the two graphs is identical,

```r
gapPermuteLambda$opt.index
[1] 29

gapERLambda$opt.index
[1] 30

GGapPermute = graph.adjacency(HugeSolutionPath$path[[gapPermuteLambda$opt.index]], mode="undirected")

GGapER = graph.adjacency(HugeSolutionPath$path[[gapERLambda$opt.index]], mode="undirected")

GapPermuteCommunities = igraph::walktrap.community(GGapPermute)

GapERCommunities = igraph::walktrap.community(GGapER)

compare(GapERCommunities, GapPermuteCommunities, method="adjusted.rand") # closer to one = closer to each other
[1] 1
```

```{r}
huge.plot(HugeSolutionPath$path[[gapPermuteLambda$opt.index]])

title("gap-com, Permuted data (pairwise correlation hard thresholding)")
```

![GapComPermuteGraph](https://user-images.githubusercontent.com/40263834/108598669-a0cd5800-7397-11eb-94e5-d63a2a6e37f9.png)

```r
huge.plot(HugeSolutionPath$path[[gapERLambda$opt.index]])

title("gap-com, ER sample (pairwise correlation hard thresholding)")
```

![GapComERGraph](https://user-images.githubusercontent.com/40263834/108598681-ae82dd80-7397-11eb-881a-5d975cd5185d.png)

The identified communities are practically identical to the ground truth clusters:

```r
TrueG = graph.adjacency(L$theta, mode = "undirected", diag = F)

TrueCommunities = igraph::walktrap.community(TrueG)

compare(GapERCommunities, TrueCommunities, method="adjusted.rand") # close to one = better
[1] 1
```

# Parallel computing

We provide an implementation which uses parallel computing. We recommend to use it when the number of variables is large:

```r
L = huge.generator(d = 1000, n = 500, graph = "cluster", g = 7)

Y = L$data

HugeSolutionPath = huge(Y, method = "ct", nlambda = nlambda)

# Without parallel:

system.time(GapLambdaER <- gap_com(HugeSolutionPath, B = 50, method = "er_sample"))
   user  system elapsed 
  87.27    1.07   88.34 

# With parallel. Use the script "gap_com_parallel.R":

library(foreach)
library(parallel)
library(doParallel)

source("gap_com_parallel.R")

registerDoParallel(cores=6)

system.time(GapLambdaERPar <- gap_com_parallel(HugeSolutionPath, B = 50, method = "er_sample"))
   user  system elapsed 
   7.17    0.60   24.42 
```

# Reference

Gap-com statistic is described in:

Kuismin and Sillanpaa (manuscript) "Gap-com: General model selection method for sparse undirected networks with clustering structure".
