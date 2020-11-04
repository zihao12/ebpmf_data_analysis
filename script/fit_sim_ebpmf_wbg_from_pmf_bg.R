## fit simulation data

library(Matrix)
library(NNLM)
library(ebpmf.alpha)
source("../code/util.R")
source("../code/init.R")
set.seed(123)

args = commandArgs(trailingOnly=TRUE)
docname = args[1]
init_name = args[2]
scale = as.integer(args[3]) ## whether to scale init
exper_version = as.integer(args[4])
maxiter = as.integer(args[5])
every = as.integer(args[6])


version = "v0.4.5"
outdir = sprintf("../output/sim/%s/exper%d", version, exper_version)
datadir = outdir
filename = sprintf("docword.%s", docname)
format = "txt"

## read data
Y = read_sim_bag_of_words(file= sprintf("%s/%s.%s",
			    datadir,filename, format))


init_file = sprintf("%s/%s_%s.Rds", datadir,docname, init_name)
init = readRDS(init_file)
init = initialize_qgl0f0w_from_l0f0LF.local(l0 = init$l0, f0 = init$f0, L = init$L, F = init$F, scale = scale)

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
	file_out = sprintf("%s/%s_ebpmf_wbg_K%d_maxiter%d_%s_scaled%d.Rds",
       outdir,docname, K, end_iter, init_name, scale)
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
