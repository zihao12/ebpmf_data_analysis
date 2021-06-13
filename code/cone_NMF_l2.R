## dim(X) = (p, n)
## init is list(B = B, W = W, R = R)
## K is the number of factors; can be null when init is provided
cone_nmf_l2 <- function(X, K, init = NULL, maxiter = 50, seed = 123){
  set.seed(seed)
  if(is.null(init)){ init <- init_cone_nmf_l2(X, K) }
  M = t(X) %*% X
  fit = cone_nmf_l2_util(M = M, B = init$B, W = init$W, R = init$R, maxiter = maxiter)
  fit[["A"]] = X %*% fit$B
  return(fit)
}


## dim(M) = (n, n); M = t(X) %*% X
## dim(B) = (n, K)
## dim(W) = (K, n); need each row has L2 norm 1
cone_nmf_l2_util <- function(M, B, W, R, maxiter = 50){
  K = ncol(B)
  loss = rep(NaN, maxiter)
  for(i in 1:maxiter){
    #browser()
    for(k in 1:K){
      Rk = R + B[,k] %o% W[k,]
      B[,k] <- threshold(Rk %*% W[k,], 1e-12)
      W[k,] <- threshold(t(Rk) %*% (M %*% B[,k]), 1e-12)
      W[k,] <- W[k,]/(sqrt(sum(W[k,]^2)))
      R <- Rk - B[,k] %o% W[k,]
    }
    loss[i] <- compute_loss_l2(M, R)
  }
  return(list(B = B, W = W, loss = loss))
}

compute_loss_l2 <- function(M, R){
  loss = sum(diag(t(R) %*% M %*% R))
  return(loss)
}

init_cone_nmf_l2 <- function(X, K){
  p = nrow(X)
  n = ncol(X)
  B = matrix(runif(n*K), nrow = n)
  W = matrix(runif(n*K), nrow = K)
  W <- W / sqrt(rowSums(W^2))
  R = diag(replicate(n, 1)) - B %*% W
  return(list(B = B, W = W, R = R))
}

threshold <- function(x, eps){
  x[x < eps] <- eps
  return(x)
}
