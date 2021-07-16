## given fitted model (a fastTopics fit for now), simulate data from multinomial
rm(list = ls())
library(fastTopics)
source("../code/simulate_multinom_util.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
modelname = args[1]
outname = args[2]
proportion = as.numeric(args[3])

# modelname = "../output/fastTopics_fit/fit_sla_fastTopics_k10.Rds"
# outname = "../data/sim_multinom/sla_simulated_fastTopics_k10_0.30.Rds"
# proportion = 0.3

model = readRDS(modelname)
if(class(model)[[1]] == "poisson_nmf_fit"){
  model = fastTopics::poisson2multinom(model)
}

sim = simulate_multinom_from_from_model(model = model, proportion = proportion)
saveRDS(sim, file = outname)
