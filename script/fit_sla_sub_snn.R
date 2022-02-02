rm(list = ls())
library(Matrix)
library(flashier)
library(magrittr)

source("../code/misc.R")

args = commandArgs(trailingOnly=TRUE)
K = as.integer(args[1])

K = 6
dat = readMM(file = "../output/fastTopics_fit/sla_small.txt")
## log-transform count data with scaling and offset 
norm.dat <- preprocess_data(dat, freq_cutoff = 1) ## won't throw away any columns of dat

## run fit with a point Laplace prior on both L and F (not imposing nonnegative constraints yet)
fit.pn <- flash.init(norm.dat, S = 1/sqrt(nrow(norm.dat)), var.type = 2) %>% flash.add.greedy(Kmax = K, prior.family = prior.point.laplace())

fit.pn <- flash.init(norm.dat, S = 1/sqrt(nrow(norm.dat)), var.type = 2
) %>% flash.init.factors(EF = fit.pn$flash.fit$EF, EF2 = fit.pn$flash.fit$EF2, prior.family = prior.point.laplace()
) %>% flash.backfit(verbose.lvl = 3, maxiter=300)

## initialize snn fit with nonnegative constraints on L, so we now have about 2*K factors
snn.init <- init.snn.LL(fit.pn)

## run snn with a point exponential prior on L, and a point Laplace prior on F, for a relatively small number of iterations 
fit.snn <- flash.init(norm.dat, S = 1/sqrt(nrow(norm.dat)), var.type = 2
) %>% flash.init.factors(EF = snn.init, prior.family = c(as.prior(ebnm::ebnm_point_exponential, sign = 1), prior.point.laplace())
) %>% flash.backfit(verbose.lvl = 3, maxiter = 25)

## refit snn with the top K factors based on pve
kset <- (length(fit.snn$pve) - rank(fit.snn$pve) < K) & (fit.snn$pve > 0)

fit.snn <- flash.init(norm.dat, S = 1/sqrt(nrow(norm.dat)), var.type = 2
) %>% flash.init.factors(EF = lapply(fit.snn$flash.fit$EF, function(x) x[, kset]), EF2 = lapply(fit.snn$flash.fit$EF2, function(x) x[, kset]),
                         prior.family = c(as.prior(ebnm::ebnm_point_exponential, sign = 1), prior.point.laplace())
) %>% flash.backfit(verbose.lvl = 3, maxiter=300)

fit.snn$sampler <- NULL
fit.snn$flash.fit <- list(EF=fit.snn$flash.fit$EF, EF2=fit.snn$flash.fit$EF2)
saveRDS(fit.snn, sprintf("../output/snn_fit/sla_small_snn_k%d.Rds", K))

