# library(pheatmap)
# library(ebpmf.alpha)
# library(ebpm)

##################################################
################ Input/Output ####################
##################################################

read_bag_of_words <- function(filename, dataname){
  if(dataname == "uci") return(read_uci_bag_of_words(filename))
  if(dataname == "sla") return(read_sla_bag_of_words(filename))
  if(dataname == "nips") return(read_nips_bag_of_words(filename))
  if(dataname == "news") return(read_news_bag_of_words(filename))
  if(dataname == "sim") return(read_sim_bag_of_words(filename))
}

## function to read data
read_uci_bag_of_words <- function(filename){
	data = read.table(file = filename, skip= 3, sep = " ",
										col.names= c("i", "j", "x"))
	out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
	return(out)
}

read_sla_bag_of_words <- function(filename){
  data = read.table(file = filename, skip= 2, sep = " ",
                    col.names= c("i", "j", "x"))
  out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
  return(out)
}

read_nips_bag_of_words <- function(filename){
  data = read.table(file = filename, skip= 2, header = FALSE, sep = " ",
                    col.names= c("i", "j", "x"))
  out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
  return(out)
}

read_news_bag_of_words <- function(filename){
  data = read.table(file = filename, skip= 2, header = FALSE, sep = " ",
                    col.names= c("i", "j", "x"))
  out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
  return(out)
}

read_sim_bag_of_words <- function(filename){
  data = read.table(file = filename, skip= 2, sep = " ",
                    col.names= c("i", "j", "x"))
  out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
  return(out)
}
##################################################
################ Initialization ##################
##################################################

# bg_prior <- function(){
# 	aL = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1e3)
#   D = length(aL)
#   g = ebpm::gammamix(pi = replicate(D, 1/D), shape = aL, scale = 1/aL)
# 	return(g)
# }
#
# initialize_qgl0f0_from_LF <- function(L, F){
# 	l0 = apply(L, 1, median)
# 	l0[l0 == 0] <- 1e-8
# 	f0 = apply(F, 1, median)
# 	f0[f0 == 0] <- 1e-8
# 	L = L/l0
# 	F = F/f0
# 	qg = ebpmf.alpha::initialize_qg_from_LF(L0 = L, F0 = F)
# 	## replace g with mixture of gamma
# 	K = ncol(L)
# 	qg$gls = replicate(K, list(bg_prior()))
# 	qg$gfs = replicate(K, list(bg_prior()))
# 	return(list(qg = qg, l0 = l0, f0 = f0))
# }








