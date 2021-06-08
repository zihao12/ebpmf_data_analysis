Y = read_uci_bag_of_words("../data/uci_BoW/docword.kos.txt")
fitted = readRDS("../output/uci_BoW/v0.3.9/kos_pmf_initLF50_K20_maxiter5000.Rds")

f_projs <- list()

for(i in 1:ncol(fitted$F)){
  f = fitted$F[,i]
  f_projs[[i]] = nnls(A = t(as.matrix(Y)), b = f)
}

saveRDS(list(f_proj = f_projs, pmf = fitted), file = "../output/kos_coneFactor.Rds")
