# # Scale each column of A so that the entries in each column sum to 1;
# # i.e., colSums(scale.cols(A)) should return a vector of ones.
# scale.cols <- function (A)
#   apply(A,2,function (x) x/sum(x))
#
# # Convert the parameters (factors & loadings) for the Poisson model to
# # the factors and loadings for the multinomial model. The return value
# # "s" gives the Poisson rates for generating the "document" sizes.
# poisson2multinom <- function (F, L) {
#   L <- t(t(L) * colSums(F))
#   s <- rowSums(L)
#   L <- L / s
#   F <- scale.cols(F)
#   return(list(F = F,L = L,s = s))
# }
library(pheatmap)


KL <- function(true,est){
  sum(ifelse(true==0,0,true * log(true/est)) + est - true)
}

JS  <- function(true,est){
  0.5*(KL(true, est) + KL(est, true))
}

RMSE <- function(true, est){
  sqrt(mean((true - est)^2))
}

log_lik <- function(X, lam){
  return(sum(dpois(x = X, lambda = lam , log = T)))
}

get_prior_summary <- function(gs, log10 = FALSE, return_matrix = FALSE){
  K = length(gs)
  phi_L = gs[[1]][["scale"]]
  idx = order(phi_L, decreasing = TRUE)
  L = length(phi_L)
  Pi = matrix(, nrow = L, ncol = K)
  for(k in 1:K){
    Pi[,k] = gs[[k]][["pi"]][idx]
  }
  rownames(Pi) = paste("phi=", round(phi_L[idx], digits = 4), sep = "")
  colnames(Pi) = paste("Topic", 1:K, sep = "")
  main = "Pi"
  if(log10){
    Pi = log10(Pi + 0.001 * min(Pi[Pi > 0]))
    main = "log10(Pi)"
  }
  pheatmap(Pi, cluster_rows=FALSE, cluster_cols=FALSE, main = main)
  if(return_matrix){return(Pi)}
}

get_multinom_from_pnmf <- function (F, L) {
  u <- colSums(F)
  F <- scale.cols(F,1/u)
  L <- scale.cols(L,u)
  s <- rowSums(L)
  L <- L / s
  return(list(F = F,L = L,s = s))
}


scale.cols <- function (A, b)
  t(t(A) * b)

## turn fitted model from pmf to a `poisson_nmf_fit` object
pmf2fastTopics_pmf <- function(fit){
  class(fit) <- c("poisson_nmf_fit","list")
  return(fit)
}
