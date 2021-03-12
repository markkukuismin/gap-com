#' Gap-com statistic
#'
#' Computes the gap-com statistic which can be applied for regularization selection of sparse undirected network estimates with clustering structure. This version is compatible with huge package (tested on v 1.2.7).
#' @param SolutionPath Solution path computed with BigQuic
#' @param B Number of reference data sets generated. Default is 50.
#' @param clustering Community detection method. Can be either "walktrap" (default), "edge_betweenness" or "fast_greedy".
#' @param w Edge weights. Default NULL.
#' @param steps The length of the random walks to perform. Default 4.
#' @param Plot Should the gap-com statistic be plotted. Default FALSE.
#' @param method Determine the number of expected clusters using (i) reference data ("permute_sample") (ii) reference graph ("er_sample", default).
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
#' library(BigQuic)
#' 
#' set.seed(46023979)
#' 
#' L = huge.generator(d=100, n=120, graph = "hub", g=5)
#' 
#' Y = L$data
#' 
#' nlambda = 50
#' 
#' Y = scale(Y)
#' 
#' S = cor(Y)
#' 
#' p = ncol(S)
#' 
#' lambda.min.ratio = 0.3
#' 
#' lambda.max = max(max(S - diag(p)), -min(S - diag(p)))
#' 
#' lambda.min = lambda.min.ratio * lambda.max
#' 
#' lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
#' 
#' BQSolutionPath = BigQuic(Y, outputFileName = "BigQuicMatrix", lambda = lambda)
#' 
#' gapLambda = gap_com(BQSolutionPath, verbose = T, Plot = T, B = 50)
#' 
#'
#' @export
#'
#' @author Markku Kuismin, Mikko J. Sillanpaa
#'
#' @references Kuismin and Sillanpaa (2020) Gap-com: General model selection method for sparse undirected networks with clustering structure

gap_com_BigQuic_parallel = function(BQSolutionPath, B = 50, clustering="walktrap", w = NULL, steps = 4, Plot=F, method="er_sample"){
  
  lambda = BQSolutionPath$lambda
  
  nlambda = length(lambda)
  
  k = rep(0, nlambda) # Nmb of clusters
  
  n = nrow(BQSolutionPath$X)
  
  p = ncol(BQSolutionPath$X)
  
  if(clustering == "walktrap") f = function(G) length(table(igraph::cluster_walktrap(G, steps = steps, weights = w)$membership))
  
  if(clustering == "edge_betweenness") f = function(G) length(table(igraph::cluster_edge_betweenness(G, weights = w)$membership))
  
  if(clustering == "fast_greedy") f = function(G) length(table(igraph::cluster_fast_greedy(G, weights = w)$membership))
  
  for(j in 1:nlambda){
    
    EdgeList = read.csv(BQSolutionPath$output_file_names[j], header=F, skip=1, sep=" ")
    
    EdgeList = as.matrix(EdgeList[ , 1:2])
    
    EdgeList = EdgeList[!is.na(EdgeList[, 1]), ]
    
    EdgeList = EdgeList[!is.na(EdgeList[, 2]), ]
    
    EdgeList = EdgeList[EdgeList[, 1] == as.integer(EdgeList[, 1]), ]
    
    EdgeList = EdgeList[EdgeList[, 2] == as.integer(EdgeList[, 2]), ]
    
    EdgeList = EdgeList[EdgeList[, 1] != 0, ]
    
    EdgeList = EdgeList[EdgeList[, 2] != 0, ]
    
    EdgeList = abs(EdgeList)
    
    G = igraph::graph.edgelist(EdgeList, directed = F)
    
    G = igraph::simplify(G)
    
    k[j] = f(G)
    
  }
  
  OutFiles = rep(NA, B*nlambda)
  
  if(method == "permute_sample"){
    
    OutFiles = foreach(i=1:B, .combine = c) %dopar% {
      
      YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
      
      BQNULLSolPath = BigQuic::BigQuic(YNULL, lambda = lambda,
                                       numthreads = BQSolutionPath$numthreads,
                                       outputFileName = paste("BQNULLSolpath", Sys.getpid(), sep=""))
      
      BQNULLSolPath$output_file_names

      
    }
    
    OutFiles = matrix(OutFiles, nrow=50, byrow=T)
    
    Expk = matrix(0, B, nlambda)
    
    for(i in 1:B){
    for(j in 1:nlambda){
      
      if(file.exists(OutFiles[i, j])){
        
        EdgeList = read.csv(OutFiles[i, j], header=F, skip=1, sep=" ")
        
        EdgeList = as.matrix(EdgeList[ , 1:2])
        
        EdgeList = EdgeList[!is.na(EdgeList[, 1]), ]
        
        EdgeList = EdgeList[!is.na(EdgeList[, 2]), ]
        
        EdgeList = EdgeList[EdgeList[, 1] == as.integer(EdgeList[, 1]), ]
        
        EdgeList = EdgeList[EdgeList[, 2] == as.integer(EdgeList[, 2]), ]
        
        EdgeList = EdgeList[EdgeList[, 1] != 0, ]
        
        EdgeList = EdgeList[EdgeList[, 2] != 0, ]
        
        EdgeList = abs(EdgeList)
        
        G = igraph::graph.edgelist(EdgeList, directed = F)
        
        G = igraph::simplify(G)
        
        Expk[i, j] = f(G)
        
      }
      
    }  
    }
      
  }
  
  if(method == "er_sample"){
    
    YNULL = apply(HugeSolPath$data, 2, function(x) x[sample(1:length(x))])
    
    BQNULLSolPath = BigQuic::BigQuic(YNULL, lambda = lambda,
                                     numthreads = BQSolutionPath$numthreads,
                                     outputFileName = "BQNULLSolpath")
    
    
    
    DummySparsity = rep(0, nlambda)
    
    for(j in 1:nlambda){
      
      EdgeList = read.table(BQNULLSolPath$output_file_names[j], header=F, skip=1)
      
      EdgeList = as.matrix(EdgeList)
      
      DummySparsity[j] = (nrow(EdgeList) - p)/(p*(p-1))
      
    }
    
    remove(list=c("YNULL"))
    
    Expk = foreach(i=1:B, .combine=c) %:% foreach(j=1:nlambda, .combine=c) %dopar% {
        
      G = igraph::erdos.renyi.game(p, DummySparsity[j] , type="gnp")
        
      f(G)

    }
    
  }
  
  if(method == "er_sample") Expk = matrix(Expk, ncol=nlambda, byrow=T)
  
  sdExpk = apply(Expk, 2, sd)
  
  skk = sdExpk * sqrt(1 + 1/B)
  
  # The next one seems odd but it must be used because
  # the errors due to the parallel computing:
  
  Expk = apply(Expk, 2, function(x) sum(x)/sum(x != 0))
  
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