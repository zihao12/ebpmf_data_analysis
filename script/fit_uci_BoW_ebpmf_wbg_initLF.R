## fit kos data

library(Matrix)
library(NNLM)
library(ebpmf.alpha)
source("../code/util.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
K = as.integer(args[2])
maxiter = as.integer(args[3])
every = as.integer(args[4])

version = "v0.4.2"
datadir = "../data/uci_BoW"
outdir = sprintf("../output/uci_BoW/%s", version)
filename = sprintf("docword.%s", docname)
format = "txt"
init_iter = 50
init_file = sprintf("%s/%s_init_nmf_K%d_iter%d.Rds",
       outdir,docname, K, init_iter)

## read data
Y = read_uci_bag_of_words(file= sprintf("%s/%s.%s",
			    datadir,filename, format))

## initialization & save file
init_tmp = readRDS(init_file)
L = init_tmp$L
F = init_tmp$F

## fit with ebpmf.alpha
T = round(maxiter/every)
init = ebpmf.alpha::initialize_qgl0f0w_from_LF(L = L, F = F)
ELBO = c()
KL = c()
RUNTIME = c()

start_time = proc.time()
for(t in 1:T){
	start_iter = 1 + (t-1)*every
	end_iter = t*every
	file_out = sprintf("%s/%s_ebpmf_wbg_initLF%d_K%d_maxiter%d.Rds",
       outdir,docname, init_iter, K, end_iter)
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
