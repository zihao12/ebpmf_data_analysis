## fit simulation data

library(Matrix)
library(NNLM)
library(ebpmf.alpha)
source("../code/util.R")
source("../code/init.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
maxiter = as.integer(args[2])
every = as.integer(args[3])


version = "v0.4.5"
exper_version = 2
outdir = sprintf("../output/sim/%s/exper%d", version, exper_version)
datadir = outdir
filename = sprintf("docword.%s", docname)
format = "txt"
init_iter = 50
init_file = sprintf("%s/truth.%s.Rds", datadir,docname)

## read data
Y = read_sim_bag_of_words(file= sprintf("%s/%s.%s",
			    datadir,filename, format))

## initialization & save file
init = readRDS(init_file)
init = initialize_qgl0f0w_from_l0f0LF.local(l0 = init$l0, f0 = init$f0, L = init$L, F = init$F)
K = length(init$w)

## fit with ebpmf.alpha
T = round(maxiter/every)
ELBO = c()
KL = c()
RUNTIME = c()

start_time = proc.time()
for(t in 1:T){
	start_iter = 1 + (t-1)*every
	end_iter = t*every
	file_out = sprintf("%s/%s_ebpmf_wbg_K%d_maxiter%d_from_truth2.Rds",
       outdir,docname, K, end_iter)
	print("##########################################")
	print(sprintf("start fitting from %d iteration", start_iter))
	fit <- ebpmf.alpha::ebpmf_wbg(X = Y, K = K,
															 pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture),
															 init = init,
															 maxiter = every, verbose  = TRUE)
	print(sprintf("finish fitting %d iteration", end_iter))
	runtime = proc.time() - start_time
	## update init, elbo
	init = list(qg = fit$qg, l0 = fit$l0, f0 = fit$f0, w = fit$w)
	ELBO = c(ELBO, fit$ELBO)
	KL = c(KL, fit$KL)
	## save file
	fit[["ELBO"]] = ELBO
	fit[["KL"]] = KL
	fit[["runtime"]] = runtime
	saveRDS(fit, file = file_out)
}
print(sessionInfo())
