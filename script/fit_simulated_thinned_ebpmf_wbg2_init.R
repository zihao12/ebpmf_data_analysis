## fit on simulated data generated from real data
## init using
rm(list = ls())
library(fastTopics)
library(ebpmf.alpha)
source("../code/misc.R")
source("../code/util.R")
source("../code/init.R")
library(Matrix)
set.seed(123)


args = commandArgs(trailingOnly=TRUE)
datafile = args[1]
initfile = args[2]
outputfile = args[3]

# datafile = "../data/sim/sla_simulated_0.10.Rds"
# initfile = "../output/fastTopics_fit/fit_sla_fastTopics_k10.Rds"
# outputfile = "../output/fastTopics_fit/fit_sla_simualted_0.10_ebpmf_wbg_k10_init_truth.Rds"

data = readRDS(datafile)
Y = as(data$Y_train, "sparseMatrix")
w_idx = data$w_idx
d_idx = data$d_idx

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

## initialize fit
init <- ebpmf.alpha::initialize_wbg2_from_LF(L = init$L, F = init$F)

## fit
print("start fitting")
start = proc.time()
fit = ebpmf.alpha::ebpmf_wbg2(X = Y, K = K, init = init, maxiter = 1000, verbose = TRUE)
fit[["runtime"]] <- proc.time() - start
saveRDS(fit, file = outputfile)

