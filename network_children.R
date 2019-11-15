# Gabriel Hoffman
# November 13, 2019
#
# With Noam.
# Given a network adjacency matrix, for each pair of nodes
# count shared children or union of children

library(clusterGeneration)
library(Matrix)
library(data.table)

# define count children function
countChildren = function( C_thresh_sp, method = c("union", "intersect") ){

	method = match.arg( method )

	resCount = lapply(1:(N-1), function(i){
		if(i %% 20 == 0 ){
			cat("\r", round(i / (N-1)*100, 1), ' %   ')
		}

		value = switch( method,
			"intersect" = C_thresh_sp[,i] %*% C_thresh_sp[,i:N], 
			"union" = colSums(C_thresh_sp[,i] + C_thresh_sp[,i:N]))

		# annotate each result with the two nodes
		data.frame( node1 = rep(i,length(value)), 
					node2 = i:N, 
					count = array(value))
		})
	resCount = data.table(do.call("rbind", resCount))
	resCount = resCount[node1 < node2,]

	resCount
}


# number of nodes
N = 5000

# generate symmatric matrix
C = genPositiveDefMat( N )$Sigma
colnames(C) = paste0('N', 1:N)
rownames(C) = paste0('N', 1:N)

# threashold matrix so that 90% of entries are 0
C_thresh = C > quantile(C[lower.tri(C)], .9)
mode(C_thresh) = "integer" # convert to integer for space
diag(C_thresh) = 0 			# a node cannot be its own parent
C_thresh[lower.tri(C_thresh)] = 0 # on the top diagonal stores data

# convert thresholded matrix to sparse matrix
C_thresh_sp = as(C_thresh, "sparseMatrix")

# notice dramatic difference in size
object.size(C)
object.size(C_thresh)
object.size(C_thresh_sp)


# get union between node i and j
i = 47
j = 49

# Compute values directly
# interestion is sum of elemente-wise product
sum(C_thresh[,i] * C_thresh[,j])

# union is sum of 
sum(C_thresh[,i] + C_thresh[,j])

# Use fast method
resCount_intersect = countChildren( C_thresh_sp, "intersect")
resCount_intersect[node1==i &node2==j,]

resCount_union = countChildren( C_thresh_sp, "union")
resCount_union[node1==i &node2==j,]




















idx = combn(N,2)

count = sapply(1:ncol(idx), function(i){
	if(i %% 10000 == 0 ){
		cat("\r", round(i / length(idx)*100, 1), ' %   ')
	}

	C_thresh[,idx[1,i]] %*% C_thresh[,idx[2,i]]

	})

res = data.frame(node1 = idx[1,],
				 node2 = idx[2,],
				 count = count)

identical(res, data.frame(resCount))






object.size(C)
object.size(C_thresh)
object.size(C_thresh_sp)

C_k = kronecker( C_thresh_sp, (C_thresh_sp), FUN='*', make.dimnames = TRUE)

 C_k[,'N1:N5']



res = matrix(1, ncol=N^2) %*% C_k
res = data.frame(value = t(as.matrix(res)))
res$node1 = array(sapply(1:N, function(x) rep(paste0('N', x), N)))
res$node2 = rep(paste0('N', 1:N), N)





