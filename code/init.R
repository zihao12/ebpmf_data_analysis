## function for initialization
library(ebpmf.alpha)

initialize_qgl0f0w_from_l0f0LF.local <- function(l0, f0, L, F){
  K = ncol(L)
  w = replicate(K, 1)
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  qg = list(qls_mean = L, qls_mean_log = log(L), kl_l = replicate(K, 0),
            qfs_mean = F, qfs_mean_log = log(F), kl_f = replicate(K, 0),
            gls = replicate(K, list(bg_prior())),
            gfs = replicate(K, list(bg_prior())))
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
  qg$gls = replicate(K, list(bg_prior()))
  qg$gfs = replicate(K, list(bg_prior()))
  return(list(qg = qg, l0 = l0, f0 = f0, w = w))
}

