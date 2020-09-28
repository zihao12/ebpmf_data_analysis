

require(igraph)

dscore = function(adj, K)
{
  # implementation of the D-SCORE (Ji and Jin, 2014). 
  # Note: When the adjacency matrix is symmetric, SCORE (2014) will be used.
  #
  # input: 
  #    adj - adjacency matrix of a connected network/graph, either
  #          symmtric (for undirected networks) or asymmtric (for directed networks)
  #    K   - the number of communities
  #
  # output:  the community lables, starting from 1, 2, ...
  # 
  # depends: igraph
  #
  # authors: Pengsheng Ji and Jiashun Jin
  
  #require(igraph)
  
  # check connectivity
  if(!is.connected(graph.adjacency(adj))) stop("The network based on the adjacency matrix is not weakly connected. Use the clusters function in the igraph package to extract the (weakly connected) giant component.")
  # make sure K > 1
  if(K < 2) stop("K must be bigger than 1.")
  
  
  # if the adj is symmetric, run SCORE
  if(sum(adj != t(adj))==0){
    x = eigen(adj, symmetric=T)$vectors
    r = x[,2:K]/x[,1]
    
    m = nrow(adj)
    threshold = log(m) 
    r[r > threshold] =  threshold
    r[r < -threshold] = -threshold
    label = kmeans(r, centers=K, nstart=20, iter=100)$cluster
    
    return(label)
  }
  else{
    # D-SCORE
    SVD = svd(adj)
    rU = as.matrix((SVD$u[,2:K]) / (SVD$u[,1]) )
    rV = as.matrix((SVD$v[,2:K]) / (SVD$v[,1]) )
    
    # indicator for the support
    a = clusters(graph.adjacency(adj %*% t(adj)))
    supU = (a$membership == which.max(a$csize))
    b = clusters(graph.adjacency(t(adj) %*% adj))
    supV = (b$membership == which.max(b$csize))
    
    # threshold the values
    m = nrow(adj)
    threshold = log(m) *1
    rU[!supU, ] = NA
    rU[rU > threshold]  = threshold
    rU[rU < -threshold] = -threshold
    
    rV[!supV, ] = NA
    rV[rV > threshold]  = threshold
    rV[rV < -threshold] = -threshold
    
    # run kmeans for the intersection of supU and supV
    UVLabel = rep(-1, m)
    UVLabel[supU & supV] = kmeans(cbind(rU[supV & supU,], rV[supU & supV, ]),
                                  centers=K, nstart=20, iter=100)$cluster
    
    UMeans=matrix(0, K, K-1)
    VMeans=UMeans
    for(i in 1:K){
      UMeans[i, ] = colMeans(as.matrix(rU[UVLabel==i, ]))
      VMeans[i, ] = colMeans(as.matrix(rV[UVLabel==i, ]))
    }
    
    # cluster the nodes in the union but not in the intersection
    for(i in 1:m){
      if(supU[i] && !supV[i]){
        d = UMeans -  matrix(1, nr=K, nc=1)  %*% rU[i,]
        d = d %*% t(d)
        d = diag(d)
        UVLabel[i] = which.min(d)              
      }
      
      if(!supU[i] && supV[i]){
        d = VMeans -  matrix(1, nr=K, nc=1)  %*% rV[i,]
        d = d %*% t(d)
        d = diag(d)
        UVLabel[i] = which.min(d)              
      }
      
    }
    
    #cluster the nodes outside the union
    while(min(UVLabel) < 0){
      tmpUVLabel = UVLabel
      for(i in 1:m){
        if(UVLabel[i] < 0){
          nLinks =   rep(-1, K)
          for(j in 1:K){
            nLinks[j] = sum(adj[i, tmpUVLabel==j]) + 
              sum(adj[tmpUVLabel==j, i])
          }
          if(max(nLinks) > 0) UVLabel[i] = which.max(nLinks)
          #cat("\n", i, ":", nLinks)
          
        }
        
      }
    }
    
    
    
    return(UVLabel)
  }
  
}


 

score = function(adj, K)
{
  # implementation of the SCORE (Jin, 2014) 
  #
  # input: 
  #    adj - symmetric adjacency matrix of a connected and undirected network/graph
  #    K   - the number of communities
  #
  # output:  the community lables, starting from 1, 2, ...
  # 
  # depends: igraph
  #
  # authors: Pengsheng Ji and Jiashun Jin
  
  #require(igraph)
  
  
  # check connectivity
  if(!is.connected(graph.adjacency(adj))) stop("The network based on the adjacency matrix is not connected. Use the clusters function in the igraph package to extract the (maximal connected) giant component.")
  # make sure K > 1
  if(K < 2) stop("K must be bigger than 1.")
  
  
  # check symmetry for adj
  if(sum(adj != t(adj))>0) stop("The adjacency matrix is not symmetric.")
  
  x = eigen(adj, symmetric=T)$vectors
  r = x[,2:K]/x[,1]
  
  m = nrow(adj)
  threshold = log(m) 
  r[r > threshold] =  threshold
  r[r < -threshold] = -threshold
  label = kmeans(r, centers=K, nstart=20, iter=100)$cluster
  
  return(label)
}




plotScree=function(adj)
{
  #adj can be symmetric or asymmetric
  m = dim(adj)
  values = eigen(adj %*% t(adj))$values 
  plot(sqrt(values[1:20]), cex=2, pch=22, type="b", lty=2, lwd=3, col="red", bg="red", xlab="", ylab="")
  
}


mat2Latex  = function(x){
  # output the matrix in latex format
  m = nrow(x)
  n = ncol(x)
  for(i in 1:m){
    for(j in 1:n){
      cat(x[i,j])
      cat("  &  ")
    }
    cat("\\\\ \n")
  }
}


##############################################


FindComponents=function(adj){
  # find all components of a network with the adjacency matrix 
  # and label them into 1, 2, ... from the largest to the smallest
  # The diag of the adj does not matter
  
  
  
  m = dim(adj)
  if(m[1] != m[2]) {
    error("adjacency matrix not  square")
    return(NA)
  }
  
  m = m[1]
  #force the diag to be 1
  diag(adj)=1
  componentLabel <- rep(0, m)
  j <- 0 # componentLabel label
  
  while(min(componentLabel)==0){
    
    j = j+1
    tmp = rowSums(adj)
    tmp = tmp * (componentLabel==0)
    jstart <- which.max(tmp)
    componentLabel[jstart] <- j
    
    # find all authors connected with jstart
    while(T){ 
      new <- colSums(adj[componentLabel == j,,drop=F ]) # drop=F keep the submatrix as a matrix
      if(sum((new>0) & (componentLabel != j)) > 0) {
        # there are new connected nodes
        componentLabel[new > 0] <- j
      }else{ 
        break
      }
      
    }
    
    #print(sum(componentLabel==j))
    
   
  }
  
  # make a new label to order them from the largest
  componentN=j 
  componentSize=rep(0,componentN)
  newLabel=rep(0,m)
  for(i in 1:componentN){
    componentSize[i]= sum(componentLabel==i)
  }
  
  tt = order(componentSize, decreasing=T) # original label
  for(i in 1:componentN){
    newLabel[componentLabel== tt[i]]  = i
  }
    
  return(newLabel)
}
  