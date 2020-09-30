## the data is already preprocessed by Peter and stored here: https://github.com/stephenslab/topics/blob/master/data/nips.RData

## I only reload and save them according to my naming habits

library(readr)
library(Matrix)
load("nips.RData")

writeMM(obj = as(counts, "sparseMatrix"),file = "docword.nips.txt")
writeLines(text = colnames(counts),con = "vocab.nips.txt")
writeLines(text = rownames(counts),con = "article.nips.txt")

