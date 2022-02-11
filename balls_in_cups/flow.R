# Gabriel Hoffman 
# Feb 10, 2022
#
# Learn flow network between proportions

# source("flow.R")

# Import packages
library(lpSolve)
library(cowplot)
library(Matrix)
library(tidyverse)
library(RColorBrewer)
library(igraph)

# scale values to be counts out of target total 
adjustScale = function(x, target=10000){

	# scale to target
	x = round(x /sum(x) * target)

	# rounding can create off by 1 error
	d = target - sum(x)

	# If target not reached, subtract difference
	# from largest category
	if( d != 0){
		i = which.max(x)
	 	x[i] = x[i] + d
	}

	x
}


# learn transition matrix give starting and ending values
getTransition = function(from, to){

	if( length(from) != length(to) ){
		stop("from and to must be same length")
	}

	# make sums the same by scaling
	from = adjustScale(from)
	to = adjustScale(to)

	len = length(from)

	# set costs, with zero costs for self transitions
	costs = matrix(1, len, len)
	diag(costs) = 0

	# set constrains to be equal
	row.signs <- rep("=", len)
	col.signs <- rep("=", len)

	# solve the linear programming transport problem
	res = lp.transport(costs, "min", row.signs, from, col.signs, to)

	# return matrix
	D = res$solution
	rownames(D) = names(from)
	colnames(D) = names(to)

	D
}


# Permutate category orderings to get average flow
# across multiple optimal solutions
getTransitionAverages = function(from, to, nperm = 1000){

	if( length(from) != length(to) ){
		stop("from and to must be same length")
	}

	n = length(from)

	if( is.null(names(from))){
		names(from) = seq(1,n)
	}

	if( is.null(names(to))){
		names(to) = seq(1,n)
	}

	# evalute paths for permutations
	res = lapply(seq(1, nperm), function(a){
		# permuate order of nodes
		i = sample(n,n)
		D = getTransition(from[i], to[i])

		# re-sort to original order
		idx = match(names(from), rownames(D))

		D[idx, idx]
	})

	# average across permutations
	D.mean = Reduce("+", res) / nperm

	# Convert mean counts to rates
	D.rate = t(apply(D.mean, 1, function(x) x/ sum(x)))

	D.rate
}

plotFromToBars = function(from, to){

	ids = names(from)

	data.frame( Category = names(from), 
		from=from/sum(from), to=to/sum(to)) %>%
		pivot_longer(cols=c("from", "to")) %>%
		ggplot(aes(Category, value, fill=name)) +
			geom_bar(stat="identity", position="dodge") +
			theme_classic() +
			theme(aspect.ratio=1) +
			ylab("Proportion") + 
			scale_y_continuous(expand=c(0, 0), limits=c(0, NA)) +
			scale_fill_brewer(palette="Set1") +
			scale_x_discrete(limits=ids)
}

plotFromToMatrix = function( D.rate ){
	
	ids = rownames(D.rate)
	
	D.rate %>%
		Matrix(sparse=TRUE) %>%
		summary %>%
		as_tibble %>%
		mutate(i = ids[i], j = ids[j]) %>%
		ggplot(aes(j,i, fill=x)) + 
			geom_tile() +
			theme_classic() +
			theme(aspect.ratio=1) +
			ylab("From") + xlab("To") +
			scale_fill_gradient2(low="white", high="navy", limits=c(0,1), name="Rate") +
			scale_x_discrete(limits=ids) +
			scale_y_discrete(limits=ids) +
			geom_text(aes(j,i, label=round(x,3)), color="red") 
}



	


plotFromToNetwork = function(D.rate){

	ids = rownames(D.rate)

	cols = brewer.pal(3, "Set1")[1:2]

	g = D.rate %>%
		Matrix(sparse=TRUE) %>%
		summary %>%	
		mutate( color = c("orange", "grey")[(i==j)+1],
				i = paste0("from_", ids[i]), 
				j = paste0("to_", ids[j]), 
				weight=x) %>%
		graph.data.frame(directed=TRUE)

	n1 = length(grep("^from_", names(V(g))))
	n2 = length(grep("^to_", names(V(g))))

	V(g)$x = c(rep(1, n1), rep(2, n2))
	V(g)$y = match(gsub("^.*_(.+)$", "\\1", names(V(g))), ids)
	V(g)$color = c(rep(cols[1], n1), rep(cols[2], n2))
	V(g)$vertex.label = gsub("^.*_(.+)$", "\\1", names(V(g)))

	plot(g, 
		edge.arrow.size = .5,
		# edge.arrow.size = E(g)$weight/4,
		edge.width=E(g)$weight*4, 
		vertex.label=V(g)$vertex.label, 
		vertex.frame.color=NA,
		vertex.label.color = "black")
}




