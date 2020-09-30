## the data is already preprocessed by Peter and stored here: https://github.com/stephenslab/topics/blob/master/data/newsgroups.RData

## I only reload and save them according to my naming habits

library(readr)
library(Matrix)
load("newsgroups.RData")

writeMM(obj = as(counts, "sparseMatrix"),file = "docword.news.txt")
writeLines(text = colnames(counts),con = "vocab.news.txt")
writeLines(text = as.vector(topics),con = "topics.news.txt")

