

#' Meta-analysis using Sidak correction
#'
#' Meta-analysis using Sidak correction evaluating the minimum 
#' p-value correcting for the number of tests.
#'
#' @param pValues array of p-values 
#' @param na.rm: a logical indicating whether missing values should be removed.
#'
#' @details
#' The Sidak method corrects a p-value based on the number of tests performed.  
#' Given a p-value p and n tests, the Sidak corrected p-value is 1-(1-min(p))^n 
#'
#' @examples
#' p = runif(10)
#'
#' meta_sidak( p )
#'
meta_sidak = function(pValues, na.rm=FALSE){

	# get the number of tests
	# if na.rm==TRUE, get the number of tests without NA values
	n = ifelse( na.rm, sum(!is.na(pValues)), length(pValues))

	# return corrected p-value
	1 - (1 - min(pValues, na.rm=na.rm))^n 
}

