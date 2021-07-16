rm(list = ls())
library(fastTopics)
source("../code/misc.R")
source("../code/util.R")
library(fastTopics)
library(Matrix)
set.seed(123)


args = commandArgs(trailingOnly=TRUE)
datafile = args[1]
outputfile = args[2]

print(outputfile)

# datafile = "../data/sim_multinom/sla_simulated_fastTopics_k4_1.Rds"
# outputfile = "../output/sim_multinom/fit_fastTopics_sla_simulated_fastTopics_k4_1.Rds"

data = readRDS(datafile)
init = fastTopics::multinom2poisson(data$truth)
init <- fastTopics::init_poisson_nmf(X = data$X, L = init$L, F = init$F)
fit <- fastTopics::fit_poisson_nmf(X = data$X, fit0 = init,
                                   numiter = 1000, method = "scd", control = list(extrapolate = TRUE))
saveRDS(fit, file = outputfile)
