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

ebpmf2fastTopics_pmf <- function(fit){
  fit = list(L = fit$l0 * fit$qg$qls_mean,
       F = fit$f0 * fit$qg$qfs_mean %*% diag(fit$w))
  return(pmf2fastTopics_pmf(fit))
}


load_model_ebpmf <- function(data_dir, data_name, method_name){
  model = readRDS(sprintf("%s/%s_%s.Rds", data_dir, data_name, method_name))
  N = length(model$ELBO)
  model.summary = list(ELBO = model$ELBO[N], KL = model$KL[N],
                       E_loglik = model$ELBO[N] + model$KL[N], runtime_iter = model$runtime[[3]]/N)
  model[["summary"]] = model.summary
  return(model)
}

load_model_pmf <- function(data_dir, data_name, method_name){
  model = readRDS(sprintf("%s/%s_%s.Rds", data_dir, data_name, method_name))
  N = length(model$log_liks)
  model.summary = list(loglik = model$log_liks[N], runtime_iter = model$runtime[[3]]/N)
  model[["summary"]] = model.summary
  return(model)
}


match_topics <- function(F1, F2, use_top = NULL){
  K = ncol(F1)
  p = nrow(F1)
  id = replicate(K, NaN)
  for(i in 1:K){
    f1 = F1[,i]
    if(is.null(use_top)){
      idx = 1:p
    }else{
      idx = sort(f1, index.return = TRUE, decreasing = TRUE)$ix[1:use_top]
    }
    dist_min = Inf
    for(j in 1:K){
      dist = sum((f1[idx] - F2[idx,j])^2)
      if(dist < dist_min){
        dist_min = dist
        matched = j
      }
    }
    id[i] <- matched
  }
  return(id)
}

match_topics_top_words <- function(F1, F2, top_words1){
  K = ncol(F1)
  id = replicate(K, NaN)
  for(i in 1:K){
    f1 = F1[,i]
    dist_min = Inf
    top_words_id = top_words1[i,]
    for(j in setdiff(1:K, id[1:i])){
      f2 = F2[,j]
      dist = sum((f1[top_words_id] - f2[top_words_id])^2)
      if(dist < dist_min){
        dist_min = dist
        matched = j
      }
    }
    id[i] = matched
  }
  return(id)
}


# multinom2poisson <- function (fit, X) {
#   F <- fit$F
#   L <- fit$L * fit$s
#   return(list(L = L, F = F))
# }

# Simulate counts from the multinomial topic model with factors F,
# loadings L and sample sizes s.
simulate_multinom_counts <- function (L, F, s) {
  n <- nrow(L)
  m <- nrow(F)
  X <- matrix(0,n,m)
  P <- tcrossprod(L,F)
  for (i in 1:n)
    X[i,] <- rmultinom(1,s[i],P[i,])
  idx <- which(colSums(X) > 0)
  return(list(L = L, F = F[idx,], X = X[, idx]))
}

## X is n by p count matrix
## prob is the probability a word is samppled
## can speed up with sparse matrix
count_subsample <- function(X, prob = 0.5){
  n = nrow(X)
  p = ncol(X)
  Y_train <- matrix(rbinom(n = n*p, size = as.matrix(X), prob = prob), nrow = n)
  Y_train <- as(Y_train, "sparseMatrix")
  Y_test <- X - Y_train
  w_idx <- which(colSums(Y_train) > 0)
  d_idx <- which(rowSums(Y_train) > 0)
  return(list(Y_train = Y_train[d_idx,w_idx], Y_test = Y_test[d_idx,w_idx],
              w_idx = w_idx, d_idx = d_idx))
}




