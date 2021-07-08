source("../code/misc.R")
source("../code/util.R")
library(Matrix)
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
prob = as.numeric(args[1])

counts = read_sla_bag_of_words("../data/SLA/docword.sla.txt")
data_sampled = count_subsample(counts, prob = prob)
saveRDS(data_sampled, sprintf("../data/sim/sla_simulated_%.2f.Rds", prob))

