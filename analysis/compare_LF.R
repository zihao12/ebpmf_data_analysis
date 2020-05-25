## script to show top words for topics
rm(list = ls())
# library(pheatmap)
# library(gridExtra)
source("../code/misc.R")
source("../code/util.R")

version = "v0.3.9"
K = 100

output_dir = sprintf("../output/uci_BoW/%s/", version)
plot_file = sprintf("../docs/kos_K%d_%s_compare_LF.pdf", K, version)
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
L_pmf = model_pmf$L; F_pmf = model_pmf$F
L_bg = model$l0 * model$qg$qls_mean; F_bg = model$f0 * model$qg$qfs_mean


## scale L, F so that colSums(F) = 1
L_pmf = L_pmf %*% diag(colSums(F_pmf))
F_pmf = F_pmf %*% diag(1/colSums(F_pmf))

L_bg = L_bg %*% diag(colSums(F_bg))
F_bg = F_bg %*% diag(1/colSums(F_bg))



plot_col = 4
plot_row = round(K/2)
plot_width = 16
plot_height = round( (plot_row/plot_col) * plot_width)

pdf(plot_file, onefile = TRUE, width = plot_width, height = plot_height)
#par(mfrow = c(plot_row, plot_col))
layout( matrix(1:(plot_col * plot_row), nrow = plot_row, byrow = TRUE))
for(k in 1:K){
  plot(L_pmf[,k], L_bg[,k], pch = 18,
       xlab = "pmf", ylab = "bg", main = sprintf("loading %d", k))
  abline(0,1,lwd=1,col="blue")
}

for(k in 1:K){
  plot(F_pmf[,k], F_bg[,k], pch = 18,
       xlab = "pmf", ylab = "bg", main = sprintf("factor %d", k))
  abline(0,1,lwd=1,col="blue")
}

dev.off()





