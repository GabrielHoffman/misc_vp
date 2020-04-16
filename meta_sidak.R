

#' Meta-analysis using Sidak correction
#'
#' Meta-analysis using Sidak correction evaluating the minimum 
#' p-value correcting for the number of tests.
#'
#' @param pValues array of p-values 
#' @param ntests number of tests to correct for.  If missing, just counts number of p-values
#' @param na.rm: a logical indicating whether missing values should be removed.
#'
#' @details
#' The Sidak method corrects a p-value based on the number of tests performed.  
#' Given a p-value p and n tests, the Sidak corrected p-value is 1-(1-min(p))^n 
#'
#' @examples
#' p = runif(10)
#'
#' # assume independent tests
#' meta_sidak( p )
#' 
#' # assume correlation structure leads to 5.4 "effective tests"
#' meta_sidak( p )
#'
meta_sidak = function(pValues, ntests, na.rm=FALSE){

	if( min(pValues) < 0 || max(pValues) > 1){
		stop("p-values must be in the range [0,1]")
	}

	if( missing(ntests) ){
		# get the number of tests
		# if na.rm==TRUE, get the number of tests without NA values
		ntests = ifelse( na.rm, sum(!is.na(pValues)), length(pValues))
	}

	# return corrected p-value
	1 - (1 - min(pValues, na.rm=na.rm))^ntests
}


# How to correct for correlation structure between features 
#----------------------------------------------------------

library(poolr) # remotes::install_github("ozancinar/poolR")
library(clusterGeneration)
n_features = 20

# get the correlation matrix between peaks
C = cov2cor(genPositiveDefMat(n_features, covMethod="unifcorrmat")$Sigma)

# Compute effective number of tests
# from the correlation matrix
ntests = meff(C, method="gao")

# get the p-values
p = runif(ncol(C))

# report minumum p-value corrected for effective number of tests
meta_sidak( p, ntests=ntests)


