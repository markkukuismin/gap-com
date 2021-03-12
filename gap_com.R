#' Gap-com statistic
#'
#' Computes the gap-com statistic which can be applied for regularization selection of sparse undirected network estimates with clustering structure. This version is compatible with huge package (tested on v 1.2.7).
#' @param SolutionPath Solution path computed with huge.
#' @param B Number of reference data sets generated. Default is 50.
#' @param clustering Community detection method. Can be either "walktrap" (default), "edge_betweenness" or "fast_greedy".
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

gap_com = function(HugeSolPath, B = 50, clustering="walktrap", w = NULL, steps = 4, Plot=F, method="er_sample", verbose=F){
  
  lambda = HugeSolPath$lambda
  
  nlambda = length(lambda)
  
  k = rep(0, nlambda) # Nmb of clusters
  
  n = nrow(HugeSolPath$data)
  
  p = ncol(HugeSolPath$data)
  
  if(clustering == "walktrap") f = function(G) cluster_walktrap(G, steps=steps, weights = w)
  
  if(clustering == "edge_betweenness") f = function(G) cluster_edge_betweenness(G, weights = w)
  
  if(clustering == "fast_greedy") f = function(G) cluster_fast_greedy(G, weights = w)
  
  for(i in 1:nlambda){
    
    A = as.matrix(HugeSolPath$path[[i]])
    
    G = graph.adjacency(A, mode = "undirected", diag=F )
    
    d = f(G)
    
    k[i] = length(table(d$membership))
    
  }
  
  Expk = matrix(0, B, nlambda)
  
  for(b in 1:B){
    
   if(method == "permute_sample"){
     
     YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
	 
     HugeBootStrapSolPath = huge(YNULL, method = HugeSolPath$method, lambda=lambda, verbose = F)
     
     for(i in 1:nlambda){
       
       A = as.matrix(HugeBootStrapSolPath$path[[i]])
       
       G = graph.adjacency(A, mode = "undirected", diag=F )
       
       d = f(G)
       
       Expk[b, i] = length(table(d$membership))
       
     }
     
   }
    
    if(method == "er_sample"){
      
      if(b == 1){
        
        YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
        
        DummySolPath = huge(YNULL, method = HugeSolPath$method, lambda=lambda, verbose = F)
        
        DummySparsity = DummySolPath$sparsity
        
        remove(list=c("YNULL", "DummySolPath"))
        
      }
      
      for(i in 1:nlambda){
        
        G = erdos.renyi.game(p, DummySparsity[i] , type="gnp")
        
        if(clustering == "walktrap") d = cluster_walktrap(G, steps=steps, weights = w)
        
        if(clustering == "edge_betweenness") d = cluster_edge_betweenness(G, weights = w)
        
        if(clustering == "fast_greedy") d = cluster_fast_greedy(G, weights = w)
        
        Expk[b, i] = length(table(d$membership))
        
      }
      
    }
    
    if(verbose){
      
      Prog = paste(c("Subsampling progress ", floor(100 * b/B), "%"), collapse = "")
      cat(Prog, "\r")
      flush.console()
      
    }
    
  }
  
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