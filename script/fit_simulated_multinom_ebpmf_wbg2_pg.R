rm(list = ls())
library(fastTopics)
source("../code/misc.R")
source("../code/util.R")
library(fastTopics)
library(ebpmf.alpha)
library(Matrix)
set.seed(123)


args = commandArgs(trailingOnly=TRUE)
datafile = args[1]
outputfile = args[2]

# datafile = "../data/sim_multinom/sla_simulated_fastTopics_k20_0.1.Rds"
# outputfile = "../output/sim_multinom/fit_ebpmf_wbg2_sla_simulated_fastTopics_k20_0.1.Rds"

data = readRDS(datafile)
K = ncol(data$truth$F)

init = fastTopics::multinom2poisson(data$truth)
init <- fastTopics::init_poisson_nmf(X = data$X, L = init$L, F = init$F)
init <- ebpmf.alpha::initialize_wbg2_from_LF(L = init$L, F = init$F)
init$qg$gfs <- replicate(K, NULL)
## fit
print("start fitting")
start = proc.time()
fit = ebpmf.alpha::ebpmf_wbg2(X = data$X, K = K, init = init,
                              pm_func = list(f = ebpm::ebpm_point1_gamma,
                                             l = ebpm::ebpm_gamma,
                                             f0 = ebpm::ebpm_gamma),
                              maxiter = 1000, verbose = TRUE)
fit[["runtime"]] <- proc.time() - start
saveRDS(fit, file = outputfile)
