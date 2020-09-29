###################################################################
#
#  This is the main code to procdue the results  in the paper by Ji 
# and Jin titled  "Analysis of Statistician Coauthorship and Citation
# Networks", to appear in the Annals of Applied Statistics.
#
####################################################################


# required packages
require(igraph) # version 0.7.1. Note there are many changes in version 1.0.0
require(bibtex)
source('functions.R') #set the working directory
require(ineq)


# read the data
authorPaperBiadj = as.matrix(read.table(file="../Data/authorPaperBiadj.txt",sep="\t", header=F))
authorList = as.matrix(read.table(file="../Data/authorList.txt",sep="\t", header=F, stringsAsFactors=F))
paperCitAdj = as.matrix(read.table("../Data/paperCitAdj.txt", header=F))
paperList = read.table("../Data/paperList.txt", sep=",", stringsAsFactors=F, header=T)


# compute adjacency matrix of the coauthorship network
coauthorAdjWeighted = authorPaperBiadj %*% t(authorPaperBiadj) 
coauthorAdj = (coauthorAdjWeighted >= 1) -0 
diag(coauthorAdj) = 0  



##########################################################################
#
# 1. Productivity, patterns and trends
#
#########################################################################

authorN = length(authorList)
paperN = dim(paperList)[1]
paperYear = paperList$year

# Fig 1 Left: papers each year
# x year; y n papers
#pdf("paperN.pdf")
par(mar=c(3,3,1,1))
plot(2003:2012, as.vector(table(paperYear)), type="o", pch=15, cex=2, lwd=4,
     main="" , ylab="", xlab="", axes=F)
axis(1,cex.axis=1.5, lwd=2, at=seq(2002, 2012, len=6))
axis(2,cex.axis=1.5, lwd=2)
#dev.off()


# Fig 1 Right: # papers per author each year
paperNPerAuthor = rep(0,10)
for( i in 2003:2012){
  nPaperT = sum((paperYear==i))
  AA= authorPaperBiadj[, paperYear==i]
  AA= rowSums(AA)
  paperNPerAuthor[i-2002] = nPaperT/sum(AA>0)
}
#pdf("nPapersPerAuthor.pdf")
par(mar=c(3,3,1,1))
plot(2003:2012, paperNPerAuthor, ylim=c(0.44,0.56 ), 
     pch=15, cex=2, type="o", lwd=4, xlab="", ylab="", axes=F)
axis(1,cex.axis=1.5, lwd=2, at=seq(2002,2012, len=6))
axis(2,cex.axis=1.5, lwd=2)
#dev.off()


# Fig 2 Left: The proportion of authors who have written more than a certain number of papers
#   x is # of papers; y is proportion. Both on log scale. 
paperNPerAuthor = rowSums(authorPaperBiadj)
Fn = ecdf(paperNPerAuthor)
x = sort(paperNPerAuthor)
y = 1 - Fn(x)
#pdf("Rplot-n-paper.pdf")
par(mar=c(3,3.5,1,.51))
plot(x[x>0 & y>0],y[x>0 & y>0], log="xy", type="b", lwd=4, pch=15, cex=2,
     xlab="", ylab="",  xlim=c(1,40), axes=F)
axis(2, lwd=2, cex.axis=1.5, at=c(.00001, 0.001, 0.01, 0.1, 1),lab=c(".00001",".001",".01",".1","1") ,las=2)
axis(1, lwd=2, cex.axis=1.5, at=c(1,2,4, 8,16,32,64))
#dev.off()


# Divide the authorship credit for each paper equally among the coauthors. For example, 
# for a two-author paper, each author is credited a half.
# The distirbution of the number of papers per author is still highly skewed.
# Lorenz curve shows that the top 10% most prolific authors
# contribute 40% of the papers, and top 20% contribute 60%.
tmp = authorPaperBiadj %*% diag(1/colSums(authorPaperBiadj)) 
paperNPerAuthor = rowSums(tmp)
ineq(paperNPerAuthor, type="Gini") # The Gini coefficient is 0.505207


# Fig 2 Right: The Lorenz curve for the number of papers each author with divided contributions.
#  x is the fraction of most prolific authors; y is the fraction of papers contributed
lc.min = Lc(paperNPerAuthor, n=rep(1,length(paperNPerAuthor)))
#pdf("Rplot-LCurve.pdf")
par(mar=c(3,2,1,1))
par(pty="s") 
plot(1-lc.min$p, 1-lc.min$L, type='l',lwd=4,
     ylim=c(0,1), xlim=c(0,1), asp=1,
     xlab="",     ylab="", yaxs='i', xaxs='i',
     axes=FALSE)
segments(0,0,1,1, lty="dashed", lwd=2)
segments(0.1, 0, 0.1, 1, lty="dashed", lwd=2  )
segments(0.2, 0, 0.2, 1, lty="dashed", lwd=2  )
segments(1, 0, 1,1, lty=1, lwd=2  )
segments(0, 1, 1,1, lty=1, lwd=2  )
axis(1,lwd=2, cex.axis=1.5)
axis(2,lwd=2, cex.axis=1.5)
#dev.off()


# Fig 3 Left: number of coauthors
#   x is number of coauthors; y is proportion of authors with coauthors more than the given number
nCoauthor = colSums(coauthorAdj)
Fn = ecdf(nCoauthor)
x = sort(unique(nCoauthor))
y = 1 - Fn(x)
#pdf("Rplot-n-coauthor.pdf")
par(pty="m") 
par(mar=c(3,3.5,2,2))
plot(x[x>0 & y>0],y[x>0 & y>0], log="xy", type="b", lwd=4, pch=15, cex=2,
     xlab="", ylab="", axes=FALSE)
axis(2, cex.axis=1.5,lwd=2, at=c(0.00001, 0.001, 0.01, 0.1,1),
     lab=c(".00001",".001",".01",".1","1"),
     las=2)
axis(1,cex.axis=1.5, lwd=2, at=c(1,4,16,64))
#dev.off()


#Fig 3 Right: The average number of coauthors for all authors who
# has published in these journals that year.
tmp = rep(0, 10)
for(i in 2003:2012){
  A = authorPaperBiadj[, paperYear==i]
  A = A %*% t(A) 
  a = diag(A)
  B = A[a>0, a>0]
  B = (B>0)-0
  diag(B)= 0
  tmp[i-2002]  = mean(rowSums(B))
}
#pdf("Rplot-avg-n-coauthor-each-year.pdf")
par(pty="m") 
par(mar=c(3,3.5,0,0.5))
plot(2003:2012,tmp, lwd=4, pch=15, cex=2, type="b",
     xlab="", ylab="", axes=FALSE)
axis(1,cex.axis=1.5, lwd=2, at=seq(2002, 2012, len=6))
axis(2, cex.axis=1.5,lwd=2, at=seq(1.5, 2.5, len=11),
     las=2)
#dev.off()


# Fig 4 Left: The Lorenz curve for the number of citation received by each paper. 
xx  = rowSums(paperCitAdj)
yy = order(xx,  decreasing=T)
paperDOI = as.character(paperList$DOI)
zz =  paperDOI[yy]
zz[1:30]
#output the list of the most cited papers
#write.table(zz, file="data/top-cited.txt", col.names=F, row.names=F)
ineq(xx, type="Gini") #0.7652821

lc.min = Lc(xx, n=rep(1,length(xx)))
par(pty="s") 
# x is Fraction of Most Cited Papers
# y Fraction of Citations
#pdf("LC-citations.pdf ")
par(mar=c(3,3,1,1))
plot(1-lc.min$p, 1-lc.min$L, type='l',lwd=4,
     ylim=c(0,1), xlim=c(0,1), #asp=1,
     xlab="",     ylab="", yaxs="i", xaxs="i", axes=F)
segments(0,0,1,1, lty="dashed")
segments(0.1, 0, 0.1, 0.99, lty="dashed"  )
segments(0.2, 0, 0.2, 0.99, lty="dashed"  )
segments(1, 0, 1, 1, lwd=2)
segments(0, 1, 1, 1, lwd=2)
axis(1, lwd=2,cex.axis=1.5)
axis(2, lwd=2,cex.axis=1.5)
#dev.off()


# Fig 4 Right: The proportions of different types of citations
#  x is year; y is the prop
# calculate the type of citation
# 1 - self, 
# 2 - coauthor
# 3 - distant
type = matrix(0, paperN, paperN)
delay =  matrix(NA, paperN, paperN)
for(i in 1:paperN){
  for(j in 1:paperN){
    if(paperCitAdj[i,j]>0){
      # find delay
      delay[i,j]= paperYear[j] - paperYear[i]
      # find type
      if(sum(authorPaperBiadj[,i] * authorPaperBiadj[,j])>0){
        type[i,j]=1
      }else if( sum(colSums(matrix(authorPaperBiadj[authorPaperBiadj[,i]>0, ], ncol=paperN)) *
                      colSums(matrix(authorPaperBiadj[authorPaperBiadj[,j]>0, ], ncol=paperN)))>0){
        type[i,j]=2
        #print(i)
        #print(j)
      }else
        type[i,j]=3
    }
  }
}
# mean delay
mean(delay[type==1]) # 2.808591
mean(delay[type==2]) # 3.356436
mean(delay[type==3]) #  3.510182
mean(delay[type>0]) #   3.302517 
# proportion
sum(type==1)/sum(type>0) # 0.2766515
sum(type==2) /sum(type>0)# 0.08825585
sum(type==3) /sum(type>0)# 0.6350926

paperCitYear = matrix(1, nrow=paperN, ncol=1) %*% t(paperYear)
prop = matrix(-1, nrow=3,ncol=5)
for(j in 1:5)
  for(i in 1:3){
    prop[i,j] = sum(type==i & (paperCitYear==2002+2*j-1 | paperCitYear==2002+2*j ))/sum(type>0 & (paperCitYear==2002+2*j-1 | paperCitYear==2002+2*j )) 
  }
#pdf("propCitation.pdf ")
par(mar=c(3,3,2,2))
plot(2002+2*(1:5), prop[1,], yaxs='i', ylim=c(0,.8), type = "o", lwd=4, col='red', ylab="", xlab='', pch="o", cex=2, axes=F)
lines(2002+2*(1:5), prop[2,], ylim=c(0,1), type = "o", lwd=4, col='green', pch=2, cex=2)
lines(2002+2*(1:5), prop[3,], ylim=c(0,1), type = "o", lwd=4, col='blue', pch=15, cex=2)
axis(1,cex.axis=2.5, lwd=3) # 1.5 and 2 for the paper
axis(2,cex.axis=2.5, lwd=3)
#legend("topleft", c("self", "coauthor", "distant"), 
#       lty=c(1,1,1), pch=c(1, 2, 15), col=c('red', 'green', 'blue'), 
#       lwd=c(4,4,4), cex=4.1,  pt.cex = 2)
# explanation: authors are getting increasingly knowledible. (tech & conferences)
#dev.off()


# reciprocate citations
citAdjWeight = authorPaperBiadj %*% paperCitAdj %*% t(authorPaperBiadj);
# force positive numbers into 1
citAdj = (citAdjWeight >0) - 0
recip = matrix(0, paperN, paperN)
for(i in 1:paperN)
  for( j in 1:paperN)
    if(paperCitAdj[i,j]>0){
      recip[i,j] = (sum(citAdj[authorPaperBiadj[,j]>0 ,  authorPaperBiadj[,i]>0])>0)
    }
sum(type==2 & recip==1) /sum(type==2) #0.7920792
sum(type==3 & recip==1) /sum(type==3) #0.2465603


# calculate fraction of the time individuals who share a
#common coauthor but have not previously collaborated
#themselves later write a paper together.
withCommonCoauthor =  coauthorAdj %*% t(coauthorAdj)
withCommonCoauthor =  (withCommonCoauthor>0)
diag(withCommonCoauthor)=0
#write.table(coauthorAdj, file="data/coauthorAdj.txt", col.names=F, row.names=F)
sum(withCommonCoauthor & !coauthorAdj) /sum(withCommonCoauthor ) #0.8119073
transitivity(graph.adjacency(coauthorAdj, mode="undirected"), type="global")
# 0.3204859
density = sum(coauthorAdj)/2/choose(authorN,2) # .00086




####################################################
#
#   centrality
#
####################################################


# Table 2:  top authors
nTop = 3
nCol = 5
authorMat = matrix(NA, nr=nTop, nc=nCol)
# # papers 
aa = rowSums(authorPaperBiadj)
bb = order(aa, decreasing=T)
authorMat[,1] = (authorList[bb[1:nTop]])
# degree (# coauthors)
aa = colSums(coauthorAdj)
bb = order(aa, decreasing=T)
authorMat[,2] = authorList[bb[1:nTop]]
# # citers (the authors who cites this person's work)
aa = rowSums(citAdj)
bb = order(aa, decreasing=T)
authorMat[,3]= (authorList[bb[1:nTop]])
# closeness
aa = closeness(graph.adjacency(coauthorAdj, mode="undirected"))
bb = order(aa, decreasing=T)
authorMat[,4]  = (authorList[bb[1:nTop]])
#betweenness
aa = betweenness(graph.adjacency(coauthorAdj, mode="undirected"))
bb = order(aa, decreasing=T)
authorMat[,5] = (authorList[bb[1:nTop]])
#mat2Latex(authorMat)
noquote(authorMat)


# Table 3:  hot papers
nTop = 5
# By # citations received (ie, papers which cite this paper)
aa = rowSums(paperCitAdj)
bb = order(aa, decreasing=T)
cat("By #Citations\n")
#bib[bb[1:nTop]]
for(i in bb[1:nTop]){
  cat(cat(authorList[authorPaperBiadj[,i]==1], sep=" and "),
      " (", paperYear[i],"). ",
  paperDOI[i], sep="")
  cat("\n")
}

# By closeness
aa = closeness(graph.adjacency(paperCitAdj, mode="directed"), mode="out")
bb = order(aa, decreasing=T)
cat("By Closeness\n")
for(i in bb[1:nTop]){
  cat(cat(authorList[authorPaperBiadj[,i]==1], sep=" and "),
      " (", paperYear[i],"). ",
      paperDOI[i], sep="")
  cat("\n")
}

# By betweenness
aa = betweenness(graph.adjacency(paperCitAdj, mode="directed"), directed=T)
bb = order(aa, decreasing=T)
cat("By Betweenness\n")
for(i in bb[1:nTop]){
  cat(cat(authorList[authorPaperBiadj[,i]==1], sep=" and "),
      " (", paperYear[i],"). ",
      paperDOI[i], sep="")
  cat("\n")
}



####################################################
#
# Coauthorship network (A)
#
###################################################


# construct the network (A)
thresh = 2  # the minimum number of coauthored papers for two authors to have an edge
# compute the adjacency matrix of the network (A)
coauthorAdjThresh = (coauthorAdjWeighted >= thresh) 
diag(coauthorAdjThresh) = 0
# find and label all components
ix = clusters(graph.adjacency(coauthorAdjThresh))
componentLabel = ix$membership
# find the labels for the largest 8 components
order(ix$csize, decreasing=T)[1:12] #output:  [1]   38  113  309  278  382  223 1036   42  175  725  996  117
# find their sizes
sort(ix$csize, decreasing=T)[1:12] #output:  [1] 236  18  14  13  10   9   9   8   8   8   8   7
#write.table(componentLabel, "data/coauthorThresh2ComponentLabel.txt", col.names=F, row.names=F)


# analyze the giant component (labeled as 38)
adj = coauthorAdjThresh[componentLabel == 38, componentLabel == 38]
author = authorList[componentLabel == 38]
# list 15 high-degree nodes
cat(author[order(rowSums(adj), decreasing=T)[1:15]], sep=", ")
# save the adjacency matrix of the giant component 
#write.table(adj, file="coauthorThresh2GiantAdj.txt", col.names=F,row.names=F)
#write.table(author, file="coauthorThresh2GiantAuthorList.txt", col.names=F,row.names=F)


# Fig 5 left: scree plot of the giant component 
plotScree(adj)


# plot the giant component
#pdf("coauthor-graph-thresh-hi-dim2.pdf", width=7, height=3.5)
par(mar=c(0,0,0,0))
set.seed(0)
vertexLabel = author
vertexLabel[apply(adj, 1, "sum")<7] =""
plot(graph.adjacency(adj, mode="undirected"), 
     vertex.color="white",
     vertex.label= vertexLabel,   vertex.label.cex=1.2, 
     vertex.size=3, vertex.label.color='red', asp=.5,
     edge.color='green', asp=.5, margin=c(0,0,0,0))
#dev.off()



# Fig 6 and 7: plot the giant component and color 
#  the nodes by the community labels

#  run SCORE for the giant component
set.seed(0)
commLabel = score(adj, 2)
# save the adjacency matrix for matlab code
write.table(adj,file="coauthorThresh2GiantAdj.txt", col.name=F,row.name=F)
# prompt to run the matlab code
readline("Please run the Matlab code MatlabCode.m, and then come back \n to press any key to continue.")
#read the community labels given by  NSC, BCPL and APL with matlab code
commLabelMatlab = data.matrix(read.table(
  "coauthorThresh2GiantCommLabelK2Matlab.txt",
  sep=",", header=F))
# put all community labels into a matrix
commLabel = cbind(commLabel, commLabelMatlab)
m = dim(commLabel)
# switch the labels 1 and 2 for SCORE and BCPL for easy plots
commLabel[,1]  = 3 - commLabel[,1]
commLabel[,2]  = 3 - commLabel[,2]
commLabel[,4]  = 3 - commLabel[,4]

# sort the columns by SCORE, APL, NSC, BCPL
commLabel = commLabel[, c(1, 4, 2,3)]
colnames(commLabel) = c("SCORE", "APL", "NSC", "BCPL")
for(i in 1:m[2]){
  #vColor =  rep("white", m[1])
  #vColor[commLabel[,i]==1] = "black"
  vShape = rep( "square", m[1])
  vShape[commLabel[,i]==1] = "circle"
  par(mar=c(0,0,0,0))
  set.seed(0)
  vertexLabel = author
  vertexLabel[apply(adj, 1, "sum")<7] =""
  plot(graph.adjacency(adj, mode="undirected"), 
       vertex.color= "red", #vColor,  ### 
       vertex.shape=vShape, 
       vertex.frame.color = "white",
       vertex.label= vertexLabel,   vertex.label.cex=1.2, 
       vertex.size=2, vertex.label.color='black', asp=.5,
       edge.color=gray(.7), asp=.5, margin=c(0,0,0,0))
  cat(colnames(commLabel)[i], "\n Press Enter to continue")
  readline()
}


# Table 4: ARI and VI for SCORE, NSC, BCPL and APL
commLabel = commLabel[, c("SCORE", "NSC", "BCPL", "APL")]
method ="adjusted.rand"
compAR = matrix(0,m[2],m[2])
for(i in 1:m[2])
  for(j in 1:m[2]){
    compAR[i,j]=compare(commLabel[,i], commLabel[,j], method=method)
  }
method ="vi"
compVI = matrix(0,m[2],m[2])
for(i in 1:m[2])
  for(j in 1:m[2]){
    compVI[i,j]=compare(commLabel[,i], commLabel[,j], method=method)
  }
for( i in 1:m[2]){
  for(j in 1:m[2]){
    if(i <= j)  cat(sprintf("%.02f", compAR[i,j]), "/", sprintf("%.02f", compVI[i,j]),sep="")
    cat(" & ")
  }
  cat(" \\\\ \n")
}


# Table 5: count the overlap of the coummunities by SCORE, NSC and APL
xx = commLabel = commLabel[, c("SCORE", "NSC", "APL")]
m = dim(xx)[2]
for(i in 1:m){
  # list all size-i subsets of 1:m
  yy = combn(1:m, i)
  n = dim(yy)[2]
  # go thru each subset
  for(j in 1:n) {
    cat(colnames(xx)[yy[,j]], sep=" $\\cap$ ")
    for(k in 1:2){
      overlap = (rowSums(as.matrix(xx[,yy[,j]]) == k) == i)
      cat(" ", sum(overlap) , " ")
      #output the name list for this cell
      #author[overlap][order(rowSums(adj)[overlap], decreasing=T)[1:20]]
    }
    cat("\n")
  }
}



# Fig 8: the stat. learning component and the dimension reduction component in one figure 
set.seed(14)  
choice = (componentLabel == 113 | componentLabel == 309)
adj = coauthorAdjThresh[choice, choice]
vertexLabel = authorList[choice]
#pdf("Rplot05-wide.pdf", width=7, height=3.5)
par(mar=c(0,2,0,1.5))
plot(graph.adjacency(adj, mode="undirected"), 
     vertex.label= vertexLabel,   vertex.label.cex=1.2, 
     vertex.size=3, vertex.label.color='red',
     edge.color='green',  asp=.5)
#dev.off()


# Table 6: the authors in "Johns Hopkins", "Duke", "Stanford", "Quantile Regression" and "Experimental Design".
# labeld as 309  278  382  223 1036  175
tempLabel = c(278, 382, 223, 1036, 175)
groupName = c("Johns Hopkins", "Duke", "Stanford", "Quantile Regression", "Experimental Design")
m = length(tempLabel)
for(i in 1:m){
  cat(groupName[i],":\n\n",sep="")
  cat(authorList[componentLabel == tempLabel[i]], sep="\n")
  cat("\n\n")
}



###########################################
#
# Coauthorship network (B)
#
##########################################


# read the adjacency matrix for the giant component  
coauthorAdj = (coauthorAdjWeighted >= 1 ) + 0
diag(coauthorAdj) = 0 
# find and label all components
ix = clusters(graph.adjacency(coauthorAdj))
# find the labels for the largest 8 components
order(ix$csize, decreasing=T)[1:10] #output:   2 113  50  87  52 130  36 114 136 141
# find their sizes
sort(ix$csize, decreasing=T)[1:10] #output: 2263   23   12   12   11   10    9    9    9    9
# find the giant component
giant = (ix$membership == which.max(ix$csize))
coauthorGiantAdj = coauthorAdj[giant, giant]
coauthorGiantAuthorList = authorList[giant]
adj = coauthorGiantAdj
author =  coauthorGiantAuthorList

# Fig 5 Central
plotScree(adj)

# run SCORE
set.seed(0)
scoreLabel = score(adj, 3)
# labels 1 for Bayes, 2 for biostat and 3 for HDDA
tmp = scoreLabel 
scoreLabel[tmp == tmp[author == "James O Berger"]] = 1
scoreLabel[tmp == tmp[author == "Joseph G Ibrahim"]] = 2
scoreLabel[tmp == tmp[author == "Peter Hall"]] = 3


# Table 7: the ARI and VI of SCORE, NSC, BCPL and APL
# read the results of the latter 3 methods from Matlab 

# save the adjacency matrix for matlab code
write.table(adj,file="coauthorGiantAdj.txt", col.name=F,row.name=F)
# prompt to run the matlab code
readline("Please run the Matlab code MatlabCode.m, and then come back \n to press any key to continue.\n")
commLabelMatlab = as.matrix(read.table(
  "coauthorGiantCommLabelK3Matlab.txt",
  sep=",", header=F))
commLabel = cbind(scoreLabel, commLabelMatlab)
colnames(commLabel)[1] = "SCORE" 
nMethods = dim(commLabel)[2]
method ="adjusted.rand"
compAR = matrix(0,nMethods,nMethods)
for(i in 1:nMethods)
  for(j in 1:nMethods){
    compAR[i,j]=compare(commLabel[,i], commLabel[,j], method=method)
  }
method ="vi"
compVI = matrix(0,nMethods,nMethods)
for(i in 1:nMethods)
  for(j in 1:nMethods){
    compVI[i,j]=compare(commLabel[,i], commLabel[,j], method=method)
  }
# print in latex format 
for( i in 1:nMethods){
  for(j in 1:nMethods){
    if(i <= j)  cat(sprintf("%.02f", compAR[i,j]), "/", sprintf("%.02f", compVI[i,j]),sep="")
    cat(" & ")
  }
  cat(" \\\\ \n")
}


# Table 8: the intersections of the communties given by SCORE, NSC and APL (columns 1,2 and 4)
# adjust labesl such that 1 for Bayes, 2 for biostat and 3 for HDDA
SCORE = commLabel[,1]
NSC = commLabel[,2]
tmp = NSC
NSC[tmp == tmp[coauthorGiantAuthorList=="James O Berger"]]  = 1
NSC[tmp == tmp[coauthorGiantAuthorList=="Joseph G Ibrahim"]]  = 2
NSC[tmp == tmp[coauthorGiantAuthorList=="Peter Hall"]]  = 3
APL = commLabel[,4]
tmp = APL
APL[tmp == tmp[coauthorGiantAuthorList=="James O Berger"]]  = 1
APL[tmp == tmp[coauthorGiantAuthorList=="Joseph G Ibrahim"]]  = 2
APL[tmp == tmp[coauthorGiantAuthorList=="Peter Hall"]]  = 3
xx = cbind(SCORE, NSC, APL)
m = dim(xx)[2]
for(i in 1:m){
  # list all size-i subsets of 1:m
  yy = combn(1:m, i)
  n = dim(yy)[2]
  # go thru each subset
  for(j in 1:n) {
    cat(colnames(xx)[yy[,j]], sep=" $\\cap$ ")
    for(k in 1:3){
      overlap = (rowSums(as.matrix(xx[,yy[,j]]) == k) == i)
      cat(" ", sum(overlap) , " ")
      #output the name list for this cell
      #author[overlap][order(rowSums(adj)[overlap], decreasing=T)[1:20]]
    }
    cat("\n")
  }
  
}


# Fig 9: plot the Bayes community with SCORE==1
ii = 1
label = author
label[colSums(adj)<9] = ""
label = label[SCORE==ii]
g = graph.adjacency(adj[SCORE==ii, SCORE==ii], mode="undirected")
#pdf("coauthor-graph-berger-3-1-vertex-size-2-wide.pdf", width=7, height=3.5)
par(mar=c(0,0,0,0))
label[label=="D Walsh"] = "Daniel Walsh"
label[label=="F Liu"] = "Fei Liu"
#for( i in 1:1){
set.seed(20) # 4 for asp=1, 5 for asp=.5
lay = layout.sphere(g)
lay[30,]=c(-.3,-.1,0)
lay[26,]=c(-.55, -.75,-.2)
plot(g, 
     vertex.label=label, vertex.label.cex=1.2, 
     vertex.size=3, vertex.label.color='red',
     edge.color='green',
     layout=lay, 
     #layout=layout.lgl(g), #20
     #layout=layout.auto(g), #10
     asp=.5
)

# cat(i)
#cat ("\n Press [enter] to continue\n")
#line <- readline()
#}
#dev.off()

# Fig 10: the Biostat community with SCORE==2
ii = 2
label = author
label[colSums(adj)<13] = ""
label = label[SCORE==ii]
label[label=="J S Marron"]="Steve Marron"
label[label=="Hao Helen Zhang"]="Helen Zhang"
label[label=="Jun S Liu"]="Jun Liu"
label[label=="Joseph G Ibrahim"] = "Joseph Ibrahim"
label[label=="Eric J Feuer"]="Eric Feuer"
label[label=="Trivellore E Raghunathan"]="Trivellore Raghunathan"
label[label=="Louise M Ryan" ]="Louise Ryan" 
g = graph.adjacency(adj[SCORE==ii, SCORE==ii], mode="undirected")
#pdf("coauthor-graph-applied-Bayesian-3-2-wide.pdf", width=7, height=3.5)
par(mar=c(0,1.3,0,0))
#for( i in 3:12){
set.seed(10)
plot(g, 
     vertex.label= label, 
     vertex.label.cex=1.2, 
     vertex.size=2, vertex.label.color='red',
     edge.color='green',
     #layout=layout.fruchterman.reingold.grid(g), #14
     #layout=layout.lgl(g), #
     layout=layout.auto(g), #10
     asp=.7
)

#cat(i)
#cat ("\n Press [enter] to continue\n")
#line <- readline()
#}
##dev.off()


# Fig 11: plot the HD community with SCORE==3
ii = 3
label = author
label[colSums(adj)<18] = ""
label = label[SCORE==ii]
label[label=="Alexandre B Tsybakov" ]="Alexandre Tsybakov" 
label[label=="Robert J Tibshirani" ]="Robert Tibshirani" 
label[label=="Ciprian M Crainiceanu"]="Ciprian Crainiceanu"
label[label=="Raymond J Carroll" ]="Raymond Carroll" 
label[label=="Wolfgang Karl Hardle"]="Wolfgang Hardle"
label[label=="Christian P Robert"  ]="Christian Robert"  
label[label=="James R Robins" ]="James Robins" 
label[label== "Lawrence D Brown"]="Lawrence Brown"
label[label=="Trevor J Hastie"  ]="Trevor Hastie" 
label[label=="Marc G Genton"] ="Marc Genton"
g = graph.adjacency(adj[SCORE==ii, SCORE==ii], mode="undirected")
#pdf("coauthor-hi-dim-3-3-complete-wide.pdf", width=7, height=5.5)
par(mar=c(0,0,0,0))
#for( i in 0:24){
set.seed(14)
lay = layout.fruchterman.reingold(g)
lay[label=="Nilanjan Chatterjee"]=c(-12,360)
lay[label=="Peter Hall"]=c(-14,450)
lay[label=="Jianqing Fan"] =c(-207, 180)
lay[label=="Lixing Zhu"]=c(306, 168)
lay[label=="Runze Li"]=c(344,-91)
lay[label=="Lawrence Brown"]=c(-10,-380)
lay[label=="Raymond Carroll"] =c(-65,-24)
lay[label=="Malay Ghosh"]  =c(313,48)
lay[label=="Marc Genton"] = c(-445,34)
lay[label=="Ciprian Crainiceanu"] = c(395,250)
plot(g, vertex.color="white",
     vertex.label= label, 
     vertex.label.cex=.9, 
     vertex.size=.7, vertex.label.color='black',
     edge.color='#AAAAFF',   #'#00AAAA',
     #layout=layout.graphopt(g), #14
     #layout=layout.auto(g), #10
     layout=lay,
     asp=.78 # 2 for .5 ; 1 for 1
)

#cat(i)
#cat ("\n Press [enter] to continue\n")
#line <- readline()
#}
#dev.off()






################################################################################
#
#  Citation network (C)   
#
########################################################################


# calculate weighted adjacency matrix where the element at 
# (i,j) is the number of citations from j to i
citAdjWeight = authorPaperBiadj %*% paperCitAdj %*% t(authorPaperBiadj)
# calculate the adjacency matrix where the element at 
# (i,j) indicates if i is ever cited by j
citAdj = (citAdjWeight >= 1 ) + 0
# ignore self-citations
diag(citAdj) = 0

# find giant component in the weak sense
ix = clusters(graph.adjacency(citAdj, mode="directed" ))
giantID = which.max(ix$csize)
citGiantAdj = citAdj[ix$membership == giantID, ix$membership == giantID]
citGiantAuthorList = authorList[ix$membership == giantID]


# Fig 5 Right
plotScree(citGiantAdj)

# run SCORE for the giant component 
set.seed(0)
tmp = dscore(citGiantAdj, 3)
dscoreLabel = tmp
# adjust the labels for convenience
dscoreLabel[tmp == tmp[citGiantAuthorList=="Yoav Benjamini"]] = 1
dscoreLabel[tmp == tmp[citGiantAuthorList=="Raymond J Carroll"]] = 2
dscoreLabel[tmp == tmp[citGiantAuthorList=="Jianqing Fan"]] = 3
# 1 multiple testing (359) including  Yoav Benjamini.
# 2 Spatial and Semi-parametric (1015) including Raymond Carroll
# 3 variable selection (1280) including Jianqing Fan



# Fig 12:  plot the eigen ratios rU and rV
K = 3
label = dscoreLabel
adj = citGiantAdj
SVD = svd(adj)
rU = (SVD$u[,2:K]) / (SVD$u[,1]) 
rV = (SVD$v[,2:K]) / (SVD$v[,1]) 
# indicator for the support
a = clusters(graph.adjacency(adj %*% t(adj)))
supU = (a$membership == which.max(a$csize))
b = clusters(graph.adjacency(t(adj) %*% adj))
supV = (b$membership == which.max(b$csize))
rU[!supU,] = NA
rV[!supV,] = NA
# thresholding
rU[rU > log(dim(citGiantAdj)[1])] = log(dim(citGiantAdj)[1])
rU[rU < -log(dim(citGiantAdj)[1])] = -log(dim(citGiantAdj)[1])
rV[rV > log(dim(citGiantAdj)[1])] = log(dim(citGiantAdj)[1])
rV[rV < -log(dim(citGiantAdj)[1])] = -log(dim(citGiantAdj)[1])
#pdf("citation-clustering-result-ru-rv.pdf", width=7, height=3.5)
par(mar=c(2,1,0,0))
par(mfrow=c(1,2), ps = 12, cex = .6, cex.main = 1)
# plot rU
pchSize = 1.5
plot(  rU[label==3 & supU & supV, ], pch='.', col='red', cex=pchSize, 
       xlab="", ylab="", xlim=c(-4,8), ylim=c(-8,2.5), type='p')
points(rU[label==2 & supU & supV, ], pch="-",col='green', cex=pchSize)
points(rU[label==1 & supU & supV, ], pch='+',col='blue', cex=pchSize)
# plot rV
plot(  rV[label==3 & supU & supV, ], pch='.', col='red', xlab="", cex=pchSize,
       ylab="", xlim=c(-4,8), ylim=c(-8,2.5), type='p')
points(rV[label==2 & supU & supV, ], pch="-", col='green', cex=pchSize)
points(rV[label==1 & supU & supV, ], pch='+', col='blue', cex=pchSize)
#dev.off()


# Fig 13: Large-Scale Multiple Testing (j=1)
# Fig 14: Variable Selection (j=3)

j = 1 # 1 for Fig 13; 3 for Fig 14
# list top authors
delta = c(23, 23, 53) # threshold for the # of citers
citGiantAuthorList[label==j  &  rowSums(citGiantAdj)> delta[j]]

aa = rowSums(citGiantAdj)
adj = citGiantAdj[label==j & aa > delta[j], label==j& aa > delta[j]]
nodeLab = citGiantAuthorList[(label==j) & aa > delta[j]]
# shorten names
nodeLab[nodeLab=="Hao Helen Zhang"] = "Helen Zhang"
nodeLab[nodeLab=="Trevor J Hastie"] = "Trevor Hastie"
nodeLab[nodeLab== "Alexandre B Tsybakov"] =  "Alexandre Tsybakov"
nodeLab[nodeLab== "Bernard W Silverman" ] =  "Bernard Silverman" 
nodeLab[nodeLab=="Jonathan E Taylor" ] = "Jonathan Taylor" 
nodeLab[nodeLab=="Joel L Horowitz"  ] = "Joel Horowitz"  
nodeLab[nodeLab=="Michael R Kosorok" ] = "Michael Kosorok" 
nodeLab[nodeLab=="R Dennis Cook" ] = "Dennis Cook" 
nodeLab[nodeLab=="Robert J Tibshirani"] = "Robert Tibshirani"
nodeLab[nodeLab=="T Tony Cai" ] = "Tony Cai" 
nodeLab[nodeLab=="Emmanuel J Candes"  ] = "Emmanuel Candes" 
nodeLab[nodeLab=="Jianhua Z Huang"]="Jianhua Huang"
nodeLab[nodeLab=="Peter J Bickel"]="Peter Bickel"
g = graph.adjacency(t(adj), mode='undirected')
seeds = c(0,0,16)  
set.seed(seeds[j]) 
#for(i in 50:125){
#  set.seed(i)
lay =  layout.fruchterman.reingold(g) 
if(j==1)  lay = layout.sphere(g)
lay[nodeLab=="Bradley Efron",] = c(-1, .5,0)
lay[nodeLab=="Aad van der Vaart"] = c(.4, .1,0)
#set.seed(21)
#lay= layout.fruchterman.reingold(g, ymax=.3, ymin=-0.01, xmin=15, xmax=17)
#lay[lay[,1] < 7, 1] = 7
#lay[lay[,2] < -6, 2] = -6
lay[nodeLab=="Marina Vannucci",] = lay[nodeLab=="Marina Vannucci",] + c(0, 9)
lay[nodeLab=="Fang Yao"] = lay[nodeLab=="Fang Yao"] + c(-3,3)
lay[nodeLab=="Bernard Silverman",] = lay[nodeLab=="Bernard Silverman",] + c(0, 3)
#lay[nodeLab=="Mohsen Pourahmadi",] = lay[nodeLab=="Mohsen Pourahmadi",] + c(0, -3)
lay[nodeLab=="Peter Hall"] = lay[nodeLab=="Peter Hall"] + c(0,-1)
lay[nodeLab=="Xihong Lin"] = lay[nodeLab=="Xihong Lin"] + c(0,1)
lay[nodeLab=="Peter Bickel"] = lay[nodeLab=="Peter Bickel"] + c(-2,0)
lay[nodeLab=="Jianqing Fan"] = lay[nodeLab=="Jianqing Fan"] + c(0,.5)
lay[nodeLab=="Nicolai Meinshausen" ] = lay[nodeLab=="Nicolai Meinshausen" ] + c(0,-1)
lay[nodeLab=="Peter Hall"] = lay[nodeLab=="Peter Hall"] + c(0,+3)
#lay[nodeLab==""] = lay[nodeLab==""] + c()
filenames = c("citation-comm-other.pdf",
              "citation-comm-variable-selection.pdf",
              "citation-comm-multiple-testing.pdf")
##pdf(filenames[j], width=7, height=7)
par(mfrow=c(1,1), mar=c(0,2,0,2))
plot(g, vertex.label=nodeLab, vertex.label.cex=1.2, 
     vertex.size=3, vertex.label.color='red',layout= lay,
     edge.color='green'
)
##dev.off()
#cat("i=", i,"\n")
#line <- readline()
#} # end of for loop




# Fig 15: the substructure of the spational/semiparametric community 
# (1015 nodes, labeled as 2,  including Carroll, Roberts)

# get the giant component of this community
adj = citGiantAdj[dscoreLabel==2,dscoreLabel==2]
author = citGiantAuthorList[dscoreLabel==2]
a  = clusters(graph.adjacency(adj))
adj  = adj[a$membership == which.max(a$csize), a$membership  == which.max(a$csize) ]
author  =  author[a$membership  == which.max(a$csize) ]

# run dscore on the giant component
kcluster = 3 # number of sub-communities
set.seed(0)
labsub  = dscore(adj, kcluster)
tmp = labsub
labsub[tmp == tmp[author=="Michael Stein"]] = 1
labsub[tmp == tmp[author=="Raymond J Carroll"]] = 2
labsub[tmp == tmp[author=="Alan Gelfan"]] = 3

# plot the high degree nodes for each sub-community
delta = c(15, 11, 12)
for( j in 1:3){
  tmp  = (labsub == j )& (rowSums(adj) > delta[j])
  #sum(tmp)
  #author[tmp]
  g = graph.adjacency(adj[tmp, tmp], mode="undirected")
  filenames  = c("citation-comm-other-sub-p-spatial-stein.pdf", 
                 "citation-comm-other-sub-np-carroll.pdf",
                 "citation-comm-other-sub-np-spatial-gelfand.pdf"
  )
  ##pdf(filenames[j], width=7, height=7)
  par(mar=c(2,3,2,3))
  seeds = c(18,12,4)
  set.seed(seeds[j]) 
  #for(i in 1:20){
  #set.seed(i)
  lay =  layout.fruchterman.reingold(g)
  tmp2=lay[,2] 
  # shorten some names
  author[author=="Hao Purdue Zhang"]="Hao Zhang"
  author[author== "Ciprian M Crainiceanu" ]= "Ciprian Crainiceanu" 
  author[author=="Jeffrey S Morris"  ]= "Jeffrey Morris" 
  author[author=="Alan H Welsh"  ]= "Alan Welsh" 
  author[author=="Trivellore E Raghunathan" ]= "Trivellore Raghunathan"
  author[author=="Raymond J Carroll"   ]= "Raymond Carroll"  
  # adjust some coordinates to avoid overlap
  if(kcluster == 3){
    #j=1 Stein
    lay[author[tmp]=="Fadoua Balabdaoui", ] =   lay[author[tmp]=="Fadoua Balabdaoui", ] + c(-5,4)
    lay[author[tmp]=="Anthony OHagan", ] =   lay[author[tmp]=="Anthony OHagan", ] + c(-5,1)
    lay[author[tmp]=="Adrian E Raftery", ] =   lay[author[tmp]=="Adrian E Raftery", ] + c(-6,1)
    
    # j=2 Carroll 
    if(j == 2) lay[tmp2 < 5, 2] = tmp2[tmp2 < 5] + 12
    lay[author[tmp]=="David Ruppert", 2] =  tmp2[author[tmp]=="David Ruppert"]+2
    lay[author[tmp]=="Nilanjan Chatterjee", 2] =  tmp2[author[tmp]=="Nilanjan Chatterjee"]+1
    lay[author[tmp]=="Jeffrey S Morris", 1] =   lay[author[tmp]=="Jeffrey S Morris", 1] - 4
    lay[author[tmp]=="D Mikis Stasinopoulos", ] =   lay[author[tmp]=="D Mikis Stasinopoulos", ] + c(-3,7)
    lay[author[tmp]=="Michael A Benjamin", ] =   lay[author[tmp]=="Michael A Benjamin", ] + c(8,-2)
    lay[author[tmp]=="Silvia Shimakura", ] =   lay[author[tmp]=="Silvia Shimakura", ] + c(8,4)
    lay[author[tmp]=="Robin Henderson", ] =   lay[author[tmp]=="Robin Henderson", ] + c(2,17)
    lay[author[tmp]=="Theo Gasser", ] =   lay[author[tmp]=="Theo Gasser", ] + c(9,17)
    lay[author[tmp]=="Robert A Rigby", ] =   lay[author[tmp]=="Robert A Rigby", ] + c(-6,9)
    lay[author[tmp]=="Brian S Caffo", ] =   lay[author[tmp]=="Brian S Caffo", ] + c(0,-2)
    
    # j=3 Gelfand 
    lay[author[tmp]=="Trivellore Raghunathan", ] =  lay[author[tmp]=="Trivellore Raghunathan", ] +c( 14, -19)
    lay[author[tmp]=="Radford M Neal", ]=lay[author[tmp]=="Radford M Neal", ]+ c(-6, 11)
    lay[author[tmp]=="Natesh Pillai", 2] =  lay[author[tmp]=="Natesh Pillai", 2] +2
    lay[author[tmp]=="David M Blei", 2] = lay[author[tmp]=="David M Blei", 2] +1
    lay[author[tmp]=="Matthew J Beal" , 2] =lay[author[tmp]=="Matthew J Beal" , 2]+.7
    lay[author[tmp]=="Robert B Gramacy", ] =  lay[author[tmp]=="Robert B Gramacy", ] +c(-6, -4)
    lay[author[tmp]=="Herbert K H Lee", ] =  lay[author[tmp]=="Herbert K H Lee", ] +c(-4, -2)
    lay[author[tmp]=="Yi Li", ] =  lay[author[tmp]=="Yi Li", ] + c(0, 4)
    lay[author[tmp]=="Alexandros Beskos", ] =  lay[author[tmp]=="Alexandros Beskos", ] + c(2, 6)
    lay[author[tmp]=="Paul Fearnhead", ] =  lay[author[tmp]=="Paul Fearnhead", ] + c(4, 5)
    
  }
  
  plot(g, 
       vertex.label=author[tmp], 
       vertex.label.cex=1.2, 
       vertex.size=3, vertex.label.color='red', layout= lay,
       edge.color='green', asp=.5)

}





# Table 9: insections of (C) and (A)
componentLabelA = componentLabel
labelC = dscoreLabel
giant  = authorList[componentLabelA==38]
machine = authorList[componentLabelA==113]
dimension = authorList[componentLabelA==309]
johns = authorList[componentLabelA==278]
duke = authorList[componentLabelA==382]
stanford = authorList[componentLabelA==223]
quant = authorList[componentLabelA==1036]
design = authorList[componentLabelA==175]
x = list(citGiantAuthorList[labelC==1],
         citGiantAuthorList[labelC==2],
         citGiantAuthorList[labelC==3],
         setdiff(authorList, citGiantAuthorList))
y = list(giant, machine, dimension, johns, duke, stanford, quant, design)
for(i in 1:4){
  for(j in 1:8){
    cat(length(intersect(x[[i]], y[[j]])), " & ")
  } 
  cat("\n")
}




# Table 10: insection of (C) and (B)

labelB = scoreLabel
# read communities of (C) by DSCORE
labelC = labelC-1
labelC[labelC==0] = labelC[labelC==0]+3
x = list(coauthorGiantAuthorList[labelB==1], 
         coauthorGiantAuthorList[labelB==2],
         coauthorGiantAuthorList[labelB==3],
         setdiff(authorList, coauthorGiantAuthorList)
)
y = list(citGiantAuthorList[labelC==1],
         citGiantAuthorList[labelC==2],
         citGiantAuthorList[labelC==3],
         setdiff(authorList, citGiantAuthorList)
)
for(i in 1:4){
  for(j in 1:4){
    overlap = intersect(x[[j]], y[[i]])
    m = length(overlap)
    cat(m, " & ")
    tmp = rep(F, length(citGiantAuthorList))
    for(k in 1:m) tmp[citGiantAuthorList==overlap[k]] = T
    degree = rowSums(citGiantAdj) [tmp]
    od = order(degree, decreasing=T)
    z = citGiantAuthorList[tmp]
    # output the names counted in each cell
    write.table(z[od], file=paste("Table-10-cell(",i,"-",j,")-",m,".txt",sep=""), col.names=F, row.names=F)
    if(i==4){
      tmp = rep(F, length(authorList))
      for(k in 1:m) tmp[authorList == overlap[k]] = T
      degree = rowSums(coauthorAdj)[tmp]
      od = order(degree, decreasing=T)
      z = authorList[tmp]
      write.table(z[od], file=paste("Table-10-cell(",i,"-",j,")-",m,".txt",sep=""), col.names=F, row.names=F)
    }
  }
  cat("\n")
  #write.table()
}





# Section 5.2.3  LNSC 
# list some high-degree nodes in each community
write.table(citGiantAdj, "citGiantAdj.txt", row.name=F,col.name=F)
readline("Please run the Matlab code MatlabCode.m, \n and come back to press any key to continue.")
labelLNSC = as.matrix(read.table("citGiantCommLabelK3Matlab.txt", sep=",",header=F))
labelLNSC = labelLNSC[,1]
aa = rowSums(citGiantAdj)
for(i in 1:3){
  bb = aa[labelLNSC == i]
  authorTmp = citGiantAuthorList[labelLNSC == i]
  authorTmp =  authorTmp[order(bb,decreasing=T)]
  cat(length(authorTmp),": ", authorTmp[1:20], sep=", ")
  cat("\n")
}












 