## fit kos data

library(Matrix)
source("../code/util.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
K = as.integer(args[2])
maxiter = as.integer(args[3])


datadir = "../data/uci_bag_of_words"
outdir = "../output/uci_bag_of_words"
filename = sprintf("docword.%s", docname)
format = "txt"
file_out = sprintf("%s/%s_ebpmf_bg_K%d_maxiter%d.Rds", 
		   outdir,docname, K, maxiter)

Y = read_uci_bag_of_words(file= sprintf("%s/%s.%s", 
			    datadir,filename, format)) 


print("fitting")
runtime <- system.time(
        fit <- ebpmf.alpha::ebpmf_bg(X = Y, K = K,
				     pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture),
				     init = NULL,
				     maxiter = maxiter, verbose  = TRUE)
)
runtime
fit[["runtime"]] = runtime
saveRDS(fit, file = file_out)

print(sessionInfo())
