## fit on simulated data generated from real data
## init using
rm(list = ls())
library(fastTopics)
source("../code/misc.R")
source("../code/util.R")
library(Matrix)
set.seed(123)


args = commandArgs(trailingOnly=TRUE)
datafile = args[1]
initfile = args[2]
outputfile = args[3]

# datafile = "../data/sim/droplet_simulated.Rds"
# initfile = "../output/fastTopics_fit/fit_droplet_fastTopics_k10.Rds"
# outputfile = "../output/fastTopics_fit/fit_droplet_simualted_fastTopics_k10_init_truth.Rds"

data = readRDS(datafile)
Y = as(data$Y_train, "sparseMatrix")
idx = data$w_idx
rm(data)

init = readRDS(initfile)
## rm the columns/rows not used in `init`
init = fastTopics::multinom2poisson(init)
init$F = init$F[w_idx,]
init$L = init$L[d_idx,]
init = fastTopics::poisson2multinom(init)
init$s <- Matrix::rowSums(Y)
init <- fastTopics::multinom2poisson(init)
init <- fastTopics::init_poisson_nmf(X = Y, L = init$L, F = init$F)
K = ncol(init$F)

fit <- fastTopics::fit_poisson_nmf(X = Y, fit0 = init, numiter = 1000, method = "scd", control = list(extrapolate = TRUE))
#fit = fastTopics::multinom2poisson(fit)

saveRDS(fit, file = outputfile)
