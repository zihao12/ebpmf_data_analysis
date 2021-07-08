source("../code/misc.R")
source("../code/util.R")
library(Matrix)
set.seed(123)

load("pbmc_68k.RData")
data_sampled = count_subsample(counts)
saveRDS(data_sampled, "../data/sim/pbmc_68k_simulated.Rds")

