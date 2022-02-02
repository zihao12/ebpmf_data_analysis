rm(list = ls())
library(Matrix)
library(fastTopics)

source("../code/misc.R")
vocab = read.csv( "../output/fastTopics_fit/sla_small_vocab.txt")[, 1]
title = readLines("../data/SLA/title.sla.txt")

# args = commandArgs(trailingOnly=TRUE)
# K = as.integer(args[1])

Ks = c(6, 10, 15, 20, 30)
for(K in Ks){
	print("####################################")
	print(sprintf("############# K = %d #############", K))
	print("####################################")
	filename = sprintf("../output/snn_fit/sla_small_snn_k%d.Rds", K)
	model = readRDS(filename)
	model = list(L = model$loadings.pm[[1]], F = model$loadings.pm[[2]], 
		d = model$loadings.scale)

	for(i in 1:K){
	  ix = sort(model$F[,i], index.return=TRUE)$ix
	  top_words = vocab[tail(ix, 8)]
	  low_words = vocab[head(ix, 8)]
	  ix = sort(model$L[,i], index.return=TRUE)$ix
	  top_docs = title[tail(ix, 8)]
	  print(sprintf("########## topic %d ##############", i))
	  print("##### top positive words:")
	  print(top_words)
	  print("##### top negative words:")
	  print(low_words)
	  print("##### top document titles:")
	  print(top_docs)
	}
}