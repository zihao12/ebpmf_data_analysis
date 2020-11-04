## function for initialization
library(ebpmf.alpha)
library(Matrix)

initialize_qgl0f0w_from_l0f0LF.local <- function(l0, f0, L, F, scale = FALSE){
  K = ncol(L)
  w = replicate(K, 1)
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  if(scale){
    sf = apply(F, 1, median)
    F = F/sf
    f0 = f0 * sf
    sl = apply(L, 1, median)
    L = L/sl
    l0 = l0 * sl
  }
  qg = list(qls_mean = L, qls_mean_log = log(L), kl_l = replicate(K, 0),
            qfs_mean = F, qfs_mean_log = log(F), kl_f = replicate(K, 0),
            gls = replicate(K, list(ebpmf.alpha::bg_prior())),
            gfs = replicate(K, list(ebpmf.alpha::bg_prior())))
  return(list(qg = qg, l0 = l0, f0 = f0, w = w))
}


initialize_qgl0f0w_from_LF.local <- function(L, F){
  K = ncol(L)
  w = replicate(K, 1)
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  l0 = apply(L, 1, median)
  #l0[l0 == 0] <- 1e-8
  f0 = apply(F, 1, median)
  #f0[f0 == 0] <- 1e-8
  L = L/l0
  F = F/f0
  qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
  ## replace g with mixture of gamma
  qg$gls = replicate(K, list(ebpmf.alpha::bg_prior()))
  qg$gfs = replicate(K, list(ebpmf.alpha::bg_prior()))
  return(list(qg = qg, l0 = l0, f0 = f0, w = w))
}


initialize_qgl0f0w_random <- function(X, K, low = 0.9, up = 1.1, seed = 123){
  set.seed(seed)
  n = nrow(X)
  p = ncol(X)
  l0 = apply(X, 1, mean)
  f0 = apply(X, 2, mean)
  l0 = l0 * (sum(X)/sum(f0)/K) / sum(l0)
  w = replicate(K, 1)
  kl_l = replicate(K, 0)
  kl_f = replicate(K, 0)
  qls_mean = matrix(runif(n*K, min = low, max = up), nrow = n, ncol = K)
  qfs_mean = matrix(runif(p*K, min = low, max = up), nrow = p, ncol = K)
  qg = list(qls_mean = qls_mean, qls_mean_log = log(qls_mean), kl_l = kl_l,
            qfs_mean = qfs_mean, qfs_mean_log = log(qfs_mean), kl_f = kl_f,
            gls = replicate(K, list(ebpmf.alpha::bg_prior())),
            gfs = replicate(K, list(ebpmf.alpha::bg_prior())))
  return(list(qg = qg, l0 = l0, f0 = f0, w = w))
}

