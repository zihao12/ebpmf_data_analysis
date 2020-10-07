library(fastTopics)
source("../code/misc.R")
source("../code/util.R")

args = commandArgs(trailingOnly=TRUE)
data_name = args[1]
K = as.integer(args[2])
maxiter = 2000
version = "v0.4.5"

dir_dict = list(sla = "SLA", nips = "NIPS", news = "News")
dir_name = dir_dict[[data_name]]
data_file = sprintf("../data/%s/docword.%s.txt", dir_name, data_name)

X = read_bag_of_words(data_file , data_name)
model = readRDS(sprintf("../output/%s/%s/%s_pmf_initLF50_K%d_maxiter%d.Rds",
												dir_name, version, data_name, K, maxiter))
diff_count <- fastTopics::diff_count_analysis(pmf2fastTopics_pmf(model), X)

saveRDS(diff_count, 
				sprintf("../output/%s/%s/%s_pmf_initLF50_K%d_maxiter%d_diffcount.Rds",
								dir_name, version, data_name, K, maxiter))
