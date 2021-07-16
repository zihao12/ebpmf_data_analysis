## functions for simulating multinomial from fitted model

normalize.cols <- function (A)
  t(t(A) / colSums(A))

simulate_multinom_counts <- function (L, F, s) {
  n <- nrow(L)
  m <- nrow(F)
  X <- matrix(0,n,m)
  P <- tcrossprod(L,F)
  for (i in 1:n)
    X[i,] <- rmultinom(1,s[i],P[i,])
  return(X)
}


simulate_multinom_from_from_model <- function(model, proportion){
  model$s <- model$s * proportion
  X = simulate_multinom_counts(L = model$L, F = model$F, s = model$s)
  w_idx <- (colSums(X) > 0)
  d_idx <- (rowSums(X) > 0)
  X <- as(X[d_idx, w_idx], "sparseMatrix")

  model$F <- normalize.cols(model$F[w_idx,])
  model$L <- model$L[d_idx,]
  model$s <- model$s[d_idx]
  return(list(X = X, truth = model))
}


