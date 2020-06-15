
ModifyBigQUIC = function(BQSolPath){
  
  lambda = BQSolPath$lambda
  
  nlambda = length(lambda)
  
  p = ncol(BQSolPath$X)
  
  path = list()
  
  for(j in 1:nlambda){
    
    EdgeList = read.csv(BQSolPath$output_file_names[j], 
                        header=F, skip=1, sep=" ") 
    
    EdgeList = as.matrix(EdgeList[ , 1:2])
    
    EdgeList = EdgeList[!is.na(EdgeList[, 1]), ]
    
    EdgeList = EdgeList[!is.na(EdgeList[, 2]), ]
    
    EdgeList = EdgeList[EdgeList[, 1] == as.integer(EdgeList[, 1]), ]
    
    EdgeList = EdgeList[EdgeList[, 2] == as.integer(EdgeList[, 2]), ]
    
    EdgeList = EdgeList[EdgeList[, 1] != 0, ]
    
    EdgeList = EdgeList[EdgeList[, 2] != 0, ]
    
    EdgeList = abs(EdgeList)
    
    A = matrix(0, p, p)
  
    A[EdgeList] = 1
    
    diag(A) = 0
    
    path[[j]] = A
      
  }
  
  obj = list(path = path, lambda = lambda)
  
  class(obj) = "huge"
  
  return(obj)
  
}