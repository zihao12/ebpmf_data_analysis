## fit kos data

library(Matrix)
library(ebpmf.alpha)
source("../code/util.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
K = as.integer(args[2])
maxiter = as.integer(args[3])
every = maxiter

version = "v0.4.5"
datadir = "../data/SLA"
outdir = sprintf("../output/SLA/%s", version)
filename = sprintf("docword.%s", docname)
format = "txt"
init_iter = 50
init_file = sprintf("%s/%s_init_nmf_K%d_iter%d.Rds",
       outdir,docname, K, init_iter)

## read data
Y = read_sla_bag_of_words(file= sprintf("%s/%s.%s",
			    datadir,filename, format))

## initialization & save file
#init = readRDS(init_file)
init = NULL

## fit with ebpmf.alpha
T = round(maxiter/every)
log_liks = c()
RUNTIME = c()

start_time = proc.time()
for(t in 1:T){
	start_iter = 1 + (t-1)*every
	end_iter = t*every
	print("##########################################")
	print(sprintf("start fitting from %d iteration", start_iter))
	fit <- ebpmf.alpha::pmf(X = Y, K = K, init = init, maxiter = every, verbose = TRUE)
	print(sprintf("finish fitting %d iteration", end_iter))
	runtime = proc.time() - start_time
	## update init, elbo
	init = list(L = fit$L, F = fit$F)
	log_liks = c(log_liks, fit$log_liks)
	## save file
	#fit[["log_liks"]] = log_liks
	#fit[["runtime"]] = runtime
	saveRDS(fit, file = init_file)
}
print(sessionInfo())
