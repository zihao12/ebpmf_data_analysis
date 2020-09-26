## script to show top words for topics
rm(list = ls())
library(pheatmap)
library(gridExtra)
source("../code/misc.R")
source("../code/util.R")

version_wbg = "v0.4.2"
version_pmf = "v0.3.9"

K = 20

output_dir_wbg = sprintf("../output/uci_BoW/%s/", version_wbg)
output_dir_pmf = sprintf("../output/uci_BoW/%s/", version_pmf)

plot_file = sprintf("../docs/kos_K%d_%s_topic_words.pdf", K, version_wbg)
data_dir = "../data/uci_BoW/"
model_wbg_initLF_name = sprintf("kos_ebpmf_wbg_initLF50_K%d_maxiter2000.Rds", K)
model_wbg_initL_name = sprintf("kos_ebpmf_wbg_initL50_K%d_maxiter2000.Rds", K)
model_pmf_name = sprintf("kos_pmf_initLF50_K%d_maxiter2000.Rds", K)
dict_name = "vocab.kos.txt"
data_name = "docword.kos.txt"

Y = read_uci_bag_of_words(file= sprintf("%s/%s",
                                        data_dir,data_name))
model_wbg1 = readRDS(sprintf("%s/%s", output_dir_wbg, model_wbg_initLF_name))
model_wbg2 = readRDS(sprintf("%s/%s", output_dir_wbg, model_wbg_initL_name))
model_pmf = readRDS(sprintf("%s/%s", output_dir_pmf, model_pmf_name))
dict = read.csv(sprintf("%s/%s", data_dir, dict_name), header = FALSE)[,1]
dict = as.vector(dict)


n = nrow(Y)
p = ncol(Y)
n_top_word = round(0.002 * p)
L_pmf = model_pmf$L
F_pmf = model_pmf$F

L_wbg1 = model_wbg1$l0 * model_wbg1$qg$qls_mean
F_wbg1 = model_wbg1$f0 * model_wbg1$qg$qfs_mean

L_wbg2 = model_wbg2$l0 * model_wbg2$qg$qls_mean
F_wbg2 = model_wbg2$f0 * model_wbg2$qg$qfs_mean

## scale L, F so that colSums(F) = 1
L_pmf = L_pmf %*% diag(colSums(F_pmf))
F_pmf = F_pmf %*% diag(1/colSums(F_pmf))

L_wbg1 = L_wbg1 %*% diag(colSums(F_wbg1))
F_wbg1 = F_wbg1 %*% diag(1/colSums(F_wbg1))

L_wbg2 = L_wbg2 %*% diag(colSums(F_wbg2))
F_wbg2 = F_wbg2 %*% diag(1/colSums(F_wbg2))


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
  F = model_wbg1$qg$qfs_mean
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("wbg_initLF topic %d", k)
  p1 = show_topic_util(F = F, idx = idx, main = main)

  F = model_wbg2$qg$qfs_mean
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("wbg_initL topic %d", k)
  p2 = show_topic_util(F = F, idx = idx, main = main)

  F = F_pmf
  idx = order(F[,k], decreasing = TRUE)[1:n_top_word]
  main = sprintf("pmf topic %d", k)
  p3 = show_topic_util(F = F, idx = idx, main = main)

  grid.arrange(grobs = list(p1, p2, p3), ncol = 2)
}

dev.off()





