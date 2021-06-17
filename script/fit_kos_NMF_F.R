rm(list = ls())
source("../code/cone_NMF_l2.R")
source("../code/util.R")
Y = read_uci_bag_of_words("../data/uci_BoW/docword.kos.txt")
init <- readRDS("../output/fitted_kos_coneNMF_F_K20.Rds")


set.seed(123)
start = proc.time()
fitted = NNLM::nnmf(A = t(as.matrix(Y)), k = 20,
                    loss = "mse", method = "scd",
                    init = list(W = init$A, H = init$W),
                    max.iter = 1000, verbose = 2)
runtime <- proc.time() - start

fitted[["runtime"]] <- runtime

saveRDS(fitted, "../output/fitted_kos_NMF_F_K20.Rds")
