## fit kos data

library(Matrix)
source("../code/util.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
K = as.integer(args[2])
maxiter = as.integer(args[3])
every = as.integer(args[4])

version = "v0.3.9"
datadir = "../data/uci_BoW"
outdir = sprintf("../output/uci_BoW/%s", version)
filename = sprintf("docword.%s", docname)
format = "txt"
file_out = sprintf("%s/%s_ebpmf_bg_K%d_maxiter%d.Rds", 
		   outdir,docname, K, maxiter)

Y = read_uci_bag_of_words(file= sprintf("%s/%s.%s", 
			    datadir,filename, format)) 



T = round(maxiter/every)
init = NULL
ELBO = c()
RUNTIME = c()

start_time = proc.time()
for(t in 1:T){
	start_iter = 1 + (t-1)*every
	end_iter = t*every
	file_out = sprintf("%s/%s_ebpmf_bg_K%d_maxiter%d.Rds",
       outdir,docname, K, end_iter)
	print("##########################################")
	print(sprintf("start fitting from %d iteration", start_iter))
	fit <- ebpmf.alpha::ebpmf_bg(X = Y, K = K,
															 pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture),
															 init = init,
															 maxiter = every, verbose  = TRUE)
	print(sprintf("finish fitting %d iteration", end_iter))
	runtime = proc.time() - start_time
	## update init, elbo
	init = list(qg = fit$qg, l0 = fit$l0, f0 = fit$f0)
	ELBO = c(ELBO, fit$ELBO)
	## save file
	fit[["ELBO"]] = ELBO
	fit[["runtime"]] = runtime
	saveRDS(fit, file = file_out)
}
print(sessionInfo())
