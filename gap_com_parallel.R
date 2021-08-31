#' Gap-com statistic
#'
#' Computes the gap-com statistic which can be applied for regularization selection of sparse undirected network estimates with clustering structure. This version is compatible with huge package (tested on v 1.2.7).
#' @param HugeSolPath Solution path computed with huge.
#' @param B Number of reference data sets generated. Default is 50.
#' @param clustering Community detection method. Can be either "walktrap" (default), "propagating_labels" or "fast_greedy".
#' @param w Edge weights. Default NULL.
#' @param steps The length of the random walks to perform. Default 4.
#' @param Plot Should the gap-com statistic be plotted. Default FALSE.
#' @param method Determine the number of expected clusters using (i) permuted data ("permute_sample") (ii) reference graph ("er_sample", default).
#' @param verbose Print the sampling progress. Default FALSE.
#' @return A list containing the following components:
#' \itemize{
#' \item opt.lambda - The largest regularization parameter value which maximizes gap-com.
#' \item opt.index - The index of the regularization parameter value which maximizes gap-com.
#' \item GapStatistic - gap-com statistic.
#' \item GapSE - Standard error of gap-com.
#' \item Expk - Expected number of network communities (clusters) under the reference distribution which is the permuted data set.
#' \item k - Estimated number of network communities.
#' \item ValidGap - Check if max(gap_com) fulfils the condition gap_com[k] >= gap_com[k+1] - GapSE[k+1]
#' \item algorithm - The community detection method used.
#' }
#' 
#' @keywords cluster gap huge network model selection sparse scio
#' @export
#' @examples
#' library(huge)
#' library(igraph)
#' library(ggplot2)
#' 
#' set.seed(46023979)
#' 
#' L = huge.generator(d=100, n=120, graph = "hub", g=5)
#' 
#' Y = L$data
#' 
#' nlambda = 50
#' 
#' HugeSolutionPath = huge(Y, method="ct", nlambda=nlambda)
#' 
#' gapLambda = gap_com(HugeSolutionPath, verbose = T, Plot = T, B = 50)
#' 
#' huge.plot(L$theta)
#' 
#' title("Ground truth")
#' 
#' huge.plot(HugeSolutionPath$path[[gapLambda$opt.index]])
#' 
#' title("gap-com")
#'
#' @export
#'
#' @author Markku Kuismin, Mikko J. Sillanpaa
#'
#' @references Kuismin and Sillanpaa (2020) Gap-com: General model selection method for sparse undirected networks with clustering structure

gap_com_parallel = function(HugeSolPath, B = 50, clustering="walktrap", w = NULL, steps = 4, Plot=F, method="er_sample"){
  
  lambda = HugeSolPath$lambda
  
  nlambda = length(lambda)
  
  k = rep(0, nlambda) # Nmb of clusters
  
  n = nrow(HugeSolPath$data)
  
  p = ncol(HugeSolPath$data)
  
  if(clustering == "walktrap") f = function(G) length(table(igraph::cluster_walktrap(G, steps = steps, weights = w)$membership))
  
  if(clustering == "propagating_labels") f = function(G) length(table(igraph::cluster_label_prop(G, weights = w)$membership))
  
  if(clustering == "fast_greedy") f = function(G) length(table(igraph::cluster_fast_greedy(G, weights = w)$membership))
  
  for(i in 1:nlambda){
    
    A = as.matrix(HugeSolPath$path[[i]])
    
    G = igraph::graph.adjacency(A, mode = "undirected", diag=F )
    
    k[i] = f(G)
    
  }
  
  Expk = rep(0, B*nlambda)
  
  d = rep(0, nlambda)
  
  if(method == "permute_sample"){
    
    Expk = foreach(i=1:B, .combine=c) %dopar% {
      
      YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
      
      HugeRefDataSolPath = huge::huge(YNULL, method = HugeSolPath$method, 
                                        lambda=lambda, verbose = F)
      
      for(j in 1:nlambda){
        
        A = as.matrix(HugeRefDataSolPath$path[[j]])
        
        G = igraph::graph.adjacency(A, mode = "undirected", diag=F )
        
        d[j] = f(G)
        
      }
      
      d
      
    }
    
  }
  
  if(method == "er_sample"){
    
    YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
    
    DummySolPath = huge(YNULL, method = HugeSolPath$method, lambda=lambda, verbose = F)
    
    DummySparsity = DummySolPath$sparsity
    
    remove(list=c("YNULL", "DummySolPath"))
    
    Expk = foreach(i=1:B, .combine=c) %:% foreach(j=1:nlambda, .combine=c) %dopar% {
        
      G = igraph::erdos.renyi.game(p, DummySparsity[j] , type="gnp")
        
      f(G)

    }
    
  }
  
  Expk = matrix(Expk, B, nlambda, byrow=T)
  
  sdExpk = apply(Expk, 2, sd)
  
  skk = sdExpk * sqrt(1 + 1/B)
  
  Expk = colMeans(Expk)
  
  Gap_lambda = Expk - k
  
  GapIndex = which.max(Gap_lambda)
  
  ind = nlambda:1
  
  ValGap = which.max(Gap_lambda) %in% ind[Gap_lambda[nlambda:2] >= Gap_lambda[(nlambda-1):1] - skk[(nlambda-1):1]]
  
  Results = list(opt.lambda = lambda[GapIndex], opt.index = GapIndex, GapStatistic = Gap_lambda, GapSE = skk, 
                 Expk = Expk, k = k, ValidGap = ValGap, algorithm = clustering)
  
  if(Plot == T){
    
    d = data.frame(x=lambda, y=Gap_lambda, SE=skk, lambda=lambda)
    
    Gp = qplot(data = d, x, y) + 
      geom_errorbar(aes(x = x, ymin = y - SE, ymax = y + SE), width = (max(lambda) - min(lambda))/20) +
      xlab(expression(lambda)) + 
      ylab(expression("Gap_lambda" %+-% "SE")) +
      ggtitle(method) +
      theme(plot.title = element_text(hjust = 0.5))
    
    print(Gp)
    
  }
  
  return(Results)
  
}