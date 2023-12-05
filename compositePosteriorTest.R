

#' @examples
#' library(muscat)
#' library(mashr)
#' library(SingleCellExperiment)
#'
#' data(example_sce)
#'
#' # create pseudobulk for each sample and cell cluster
#' pb <- aggregateToPseudoBulk(example_sce[1:100, ],
#'   assay = "counts",
#'   cluster_id = "cluster_id",
#'   sample_id = "sample_id",
#'   verbose = FALSE
#' )
#'
#' # voom-style normalization
#' res.proc <- processAssays(pb, ~group_id)
#'
#' # Differential expression analysis within each assay,
#' # evaluated on the voom normalized data
#' res.dl <- dreamlet(res.proc, ~group_id)
#'
#' # run MASH model
#' # This can take 10s of minutes on real data
#' # This small datasets should take ~30s
#' res_mash <- run_mash(res.dl, "group_idstim")
#' 
#' # Composite test based on posterior probabilities
#' # to identify effect present in *at least 1* monocyte type
#' # and *NO* T-cell type.
#' include = c("CD14+ Monocytes", "FCGR3A+ Monocytes")
#' exclude = c('CD4 T cells', 'CD8 T cells')
#' 
#' prob = compositePosteriorTest(res_mash, include, exclude)
#' 
#' # examine the lFSR for top gene
#' get_lfsr(x$model)[which.max(prob),,drop=FALSE]
#
#' @export
compositePosteriorTest = function( x, include, exclude, test = c("at least 1", "all")){
	
	test = match.arg(test)

	stopifnot(is(x, "dreamlet_mash_result"))

	# get probability from lFSR
	prob = 1 - get_lfsr(x$model)
	
	.compositePosteriorTest( prob, include, exclude, test)
}

# Given matrix of posterior probabilities with genes on rows and columns as conditions, compute composite probability from include vs exclude set.
.compositePosteriorTest = function( prob, include, exclude, test = c("at least 1", "all")){
	
	test = match.arg(test)

	# if probability is NA, set it to zero
	prob[is.na(prob)] = 0

	stopifnot( all(include %in% colnames(prob)) )
	stopifnot( all(exclude %in% colnames(prob)) )

	# get probability of being correct sign
	prob = 1 - get_lfsr(x$model)

	# probability that *NO* cell types have a non-zero effect
	prob_excl = apply(1 - prob[,exclude,drop=FALSE], 1, prod, na.rm=TRUE)

	if( test == "at least 1"){
		# Probability at least 1, (i.e. probability that not none)
		prob_incl = 1 - apply(1 - prob[,include,drop=FALSE], 1, prod, na.rm=TRUE)
	}else if( test == "all"){
		# for each gene
		# probability that *ALL* cell types have a non-zero effect
		prob_incl = apply(prob[,include], 1, prod, na.rm=TRUE)
	}

	prob_incl * prob_excl
}




