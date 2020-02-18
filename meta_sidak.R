

#' Meta-analysis using Sidak correction
#'
#' Meta-analysis using Sidak correction evaluating the minimum 
#' p-value correcting for the number of tests.
#'
#' @param pValues array of p-values 
#'
#' @details
#' The Sidak method corrects a p-value based on the number of tests performed.  
#' Given a p-value p and n tests, the Sidak corrected p-value is 1-(1-min(p))^n 
#'
meta_sidak = function(pValues){

	n = length(pValues)
	p.corrected = 1 - (1 - min(pValues))^n 
	p.corrected
}



