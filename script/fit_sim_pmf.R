## fit simulation data

library(Matrix)
library(NNLM)
library(ebpmf.alpha)
source("../code/util.R")
source("../code/init.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
K = as.integer(args[2])
maxiter = as.integer(args[3])
every = as.integer(args[4])
exper_version = as.integer(args[5])


version = "v0.4.5"
outdir = sprintf("../output/sim/%s/exper%d", version, exper_version)
datadir = outdir
filename = sprintf("docword.%s", docname)
format = "txt"

## read data
Y = read_sim_bag_of_words(file= sprintf("%s/%s.%s",
			    datadir,filename, format))

## initialization & save file
init = NULL
## fit with ebpmf.alpha
T = round(maxiter/every)
logliks <- c()
RUNTIME = c()

start_time = proc.time()
for(t in 1:T){
	start_iter = 1 + (t-1)*every
	end_iter = t*every
	file_out = sprintf("%s/%s_pmf_K%d_maxiter%d.Rds",
       outdir,docname, K, end_iter)
	print("##########################################")
	print(sprintf("start fitting from %d iteration", start_iter))
	fit <- ebpmf.alpha::pmf(X = Y, K = K,
															 init = init,
															 maxiter = every, verbose  = TRUE)
	print(sprintf("finish fitting %d iteration", end_iter))
	runtime = proc.time() - start_time
	## update init, elbo
	init = list(L = fit$L, F = fit$F)
	logliks = c(logliks, fit$log_liks)
	## save file
	fit[["log_liks"]] = logliks
	fit[["runtime"]] = runtime
	saveRDS(fit, file = file_out)
}
print(sessionInfo())
