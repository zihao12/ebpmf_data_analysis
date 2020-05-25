## script to show top words for topics
rm(list = ls())
library(pheatmap)
library(gridExtra)
source("../code/misc.R")
source("../code/util.R")

version = "v0.3.9"
K = 100

output_dir = sprintf("../output/uci_BoW/%s/", version)
plot_file = sprintf("../docs/kos_K%d_%s_topic_words.pdf", K, version)
data_dir = "../data/uci_BoW/"
model_name = sprintf("kos_ebpmf_bg_initLF50_K%d_maxiter2000.Rds", K)
model_pmf_name = sprintf("kos_pmf_initLF50_K%d_maxiter2000.Rds", K)
dict_name = "vocab.kos.txt"
data_name = "docword.kos.txt"

Y = read_uci_bag_of_words(file= sprintf("%s/%s",
                                        data_dir,data_name))
model = readRDS(sprintf("%s/%s", output_dir, model_name))
model_pmf = readRDS(sprintf("%s/%s", output_dir, model_pmf_name))
dict = read.csv(sprintf("%s/%s", data_dir, dict_name), header = FALSE)[,1]
dict = as.vector(dict)


n = nrow(Y)
p = ncol(Y)
n_top_word = round(0.002 * p)
L_pmf = model_pmf$L; F_pmf = model_pmf$F
L_bg = model$l0 * model$qg$qls_mean; F_bg = model$f0 * model$qg$qfs_mean
lf = poisson2multinom(L=L_bg,F=F_bg)
lf_pmf = poisson2multinom(L = L_pmf,F = F_pmf)

## scale L, F so that colSums(F) = 1
L_pmf = L_pmf %*% diag(colSums(F_pmf))
F_pmf = F_pmf %*% diag(1/colSums(F_pmf))

L_bg = L_bg %*% diag(colSums(F_bg))
F_bg = F_bg %*% diag(1/colSums(F_bg))


show_topic_util <- function(F, idx, main){
  F_sub = F[idx,]
  rownames(F_sub) = dict[idx]
  #colnames(F_sub) = paste("Topic", 1:K, sep = "")
  colnames(F_sub) = NULL
  p = pheatmap(F_sub,
           cluster_rows=FALSE, cluster_cols=FALSE,
           silent = TRUE,
           main = main)[[4]]
  return(p)
}

pdf(plot_file, onefile = TRUE)
for(k in 1:K){
  #k = 1
  ## show top words(bar_f_jk) in bar_f_jk
  F = model$qg$qfs_mean
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("bg topic %d", k)
  p1 = show_topic_util(F = F, idx = idx, main = main)
  ## show top words(bar_f_jk) in f_jk
  F = F_pmf
  main = sprintf("pmf topic %d", k)
  p2 = show_topic_util(F = F, idx = idx, main = main)
  ## show top words(f_j0 * bar_f_jk) in bar_f_jk
  F = F_bg
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("bg topic %d", k)
  p3 = show_topic_util(F = F, idx = idx, main = main)
  ## show top words (pmf) in f_jk
  F = F_pmf
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("pmf topic %d", k)
  p4 = show_topic_util(F = F, idx = idx, main = main)
  ## show top words (pmf) in f_jk
  F = lf$F
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("bg topic %d", k)
  p5 = show_topic_util(F = F, idx = idx, main = main)
  ## show top words (pmf) in f_jk
  F = lf_pmf$F
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("pmf topic %d", k)
  p6 = show_topic_util(F = F, idx = idx, main = main)

  grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6), ncol = 2)
}

dev.off()





