## function to read data

read_uci_bag_of_words <- function(filename){
	data = read.table(file = filename, skip= 3,, sep = " ",
										col.names= c("i", "j", "x"))
	out = Matrix::sparseMatrix(i = data$i, j = data$j, x = data$x)
	return(out)
}
