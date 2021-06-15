library(mvtnorm)

normalize.cols <- function (A)
  t(t(A) / colSums(A))

simulate_factors <- function (m, k) {
  n_top = 50
  F <- matrix(runif(n = m*k),m,k)
  for(i in 1:k){
    idx = sample(1:m, n_top, replace = FALSE)
    F[idx, i] <- 10 * runif(n_top)
  }
  return(normalize.cols(F))
}




# Randomly generate an n x k loadings matrix (i.e., mixture
# proportions matrix) for a multinomial topic model with k topics. The
# loadings are simulated from the correlated topic model with mean
# zero and k x k covariance matrix S.
simulate_loadings <- function (n, k, S) {
  L <- matrix(0,n,k)
  for (i in 1:n) {
    u <- rmvnorm(1,sigma = S)
    u <- u - max(u)
    L[i,] <- exp(u)/sum(exp(u))
  }
  return(L)
}

# Randomly generate samples sizes for the multinomial data simulation.
simulate_sizes <- function (n, doc_len)
  rpois(n = n, lambda = doc_len)

# Simulate counts from the multinomial topic model with factors F,
# loadings L and sample sizes s.
simulate_multinom_counts <- function (L, F, s) {
  n <- nrow(L)
  m <- nrow(F)
  X <- matrix(0,n,m)
  P <- tcrossprod(L,F)
  for (i in 1:n)
    X[i,] <- rmultinom(1,s[i],P[i,])
  return(X)
}

## L, F are truth (in multinomial model)
## model is a list containing fitted L, F (also in multinomial model)
compare_truth_fitted <- function(model, idx, L, F, n_sample){
  k = ncol(F)
  for(i in 1:k){
    f_fitted = model$F[,idx[i]]
    f_true = F[,i]
    mse = mean((f_fitted - f_true)^2)
    plot(f_fitted, f_true, main = sprintf("F topic %d: mse %f", i, mse))
    abline(a = 0, b = 1, col = "blue")
    l_fitted = model$L[1:n_sample,idx[i]]
    l_true = L[1:n_sample,i]
    mse = mean((l_fitted - l_true)^2)
    plot(l_fitted, l_true, main = sprintf("L topic %d: mse %f", i, mse))
    abline(a = 0, b = 1, col = "blue")
  }
}






