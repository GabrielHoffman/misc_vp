# Gabriel Hoffman
# December 13, 2020
#
# Compute the residual degrees of freeom for a linear mixed model

# trace only gives correct rdf if evals are 1 and 0
# trace is given in the book with H is idempotent.
# I = diag(1, n)
# tr((I-H) %*% (I-H))
# When H is not idempotent, this actually gives the degrees of freedom
# of the best chi-square fit to the distribution
# if x ~ sum(lambda * rchisq(n)), then best chi-square approximation
# is a 1 momoent match to the mean of sum(lambda)
# When no strinkage, rdf.trace is exact.
# The approximation improves as sum(lambda) increases
# For smaller rdf, the chi-square is shifted slighly left
rdf.trace = function(fit){
  tr = function(A) sum(diag(A))
   # get hat matrix from linear mixed model
  H = lme4:::hatvalues.merMod(fit, fullHatMatrix=TRUE)   

  # number of samples
  n = nrow(H) 

  # I = diag(1, n)
  # tr((I-H) %*% (I-H))
  n - 2*tr(H) + sum(H*H)
}

#' Fast approximate residual degrees of freedom 
#'
#' Defining \eqn{H = A^TA + B^TB} where \eqn{A} and \eqn{B} are low rank, compute 
#' \eqn{n - 2tr(H) + tr(HH)} in \eqn{O(np^2)} instead of \eqn{O(n^2p^2)}.
#' 
#' @param A a \code{matrix} or \code{sparseMatrix}
#' @param B a \code{matrix} or \code{sparseMatrix}
#' 
#' @examples
#' 
#' # sample size
#' n = 500
#' # number of covariates
#' n_a = 20
#' n_b = 20
#' 
#' # Simulate low rank matrices
#' A = matrix(rnorm(n_a*n), ncol=n)
#' B = matrix(rnorm(n_b*n), ncol=n)
#' 
#' # Evaluate RDF
#' rdf_from_matrices(A, B)
#' 
#' @importFrom Matrix
#' @importFrom sparsesvd sparsesvd
#' @seealso rdf.merMod
#'
rdf_from_matrices = function(A,B){

	# desired result, but this simple code is O(n^2)
	# H = crossprod(A) + crossprod(B)
	# n = nrow(H)
	# rdf = n - 2*tr(H) + sum(H*H)

	# Compute SVD of A.
	# if A is a sparseMatrix, sparsesvd() is substantially faster
	if( is(A, "sparseMatrix") ){
		dcmp_A = sparsesvd(A, rank=nrow(A))
	}else{
		dcmp_A = svd(A, nu=0)
	}

	# B
	if( is(B, "sparseMatrix") ){
		dcmp_B = sparsesvd(B, rank=nrow(B))
	}else{
		dcmp_B = svd(B, nu=0, nv=0)
	}

	# system.time(replicate(1000, svd(A, nu=0)))
	# system.time(replicate(1000, sparsesvd(A, rank=nrow(A))))

	# system.time(replicate(100, svd(B, nu=0)))
	# system.time(replicate(100, sparsesvd(B, rank=nrow(B))))

	# U = dcmp_A$u
	V = dcmp_A$v
	s_a = dcmp_A$d
	s_b = dcmp_B$d
	# Q = crossprod(U, B) %*% dcmp_A$v # change order of this?

	# Q = crossprod(U, B %*% dcmp_A$v) # change order of this?
	# sum(s_a^4) + 2*tr( diag(s_a^2) %*% crossprod(Q)) + sum(s_b^4) 
	# sum(s_a^4) + 2*sum( diag(s_a^2) * crossprod(B %*% V)) + sum(s_b^4) 
	tr_H_H = sum(s_a^4) + 2*tr( (s_a^2) * crossprod(B %*% V)) + sum(s_b^4) 

	n = ncol(A)

	n - 2 * (sum(s_a^2) + sum(s_b^2)) + tr_H_H
}


#' Approximate residual degrees of freedom
#'
#' Compute the approximate residual degrees of freedom from a linear mixed model.
#'
#' @param model An object of class \code{merMod}
#'
#' @description 
#' For a linear model with \eqn{n} samples and \eqn{p} covariates, \eqn{RSS/sigma^2 \sim \chi^2_{\nu}} where \eqn{\nu = n-p} is the residual degrees of freedom.  In the case of a linear mixed model, the distribution is no longer exactly a chi-square distribution, but can be approximated with a chi-square distribution. 
#'
#' Given the hat matrix, \code{H}, that maps between observed and fitted responses, the approximate residual degrees of freedom is \eqn{\nu = tr((I-H)(I-H))}.  For a linear model, this simplifies to the well known form \eqn{\nu = n - p}. In the more general case, such as a linear mixed model, the original form simplifies only to \eqn{n - 2tr(H) + tr(HH)} and is an approximation rather than being exact.  The third term here is quadratic time in the number of samples, \eqn{n} and can be computationally expensive to evalaute for larger datasets.  Here we use a linear time algorithm that takes advantage of the fact that \eqn{H} is low rank.
#'
#' \eqn{H} is computed as \eqn{A^TA + B^TB} for \code{A=CL} and \code{B=CR} defined in the code.  Since \eqn{A} and \eqn{B} are low rank, there is no need to compute \eqn{H} directly.  Instead, the terms \eqn{tr(H)} and \eqn{tr(HH)} can be computed using only the singular value decomposition (SVD) of \eqn{A} and \eqn{B} in \eqn{O(np^2)}. Moreover, the SVD computation can take advantage of the fact that \eqn{A} is a \code{sparseMatrix}.  
#'
#' @examples
#' # Fit linear mixed model
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # Evaluate the approximate residual degrees of freedom
#' rdf.merMod(fit)
#'
#' @import Matrix lme4
#' @seealso rdf_from_matrices
#' @export
rdf.merMod = function(model) {
	# code adapted from lme4:::hatvalues.merMod
    if (isGLMM(model)) 
        warning("the hat matrix may not make sense for GLMMs")

    sqrtW <- Diagonal(x = sqrt(weights(model, type = "prior")))

    rdf <- with(getME(model, c("L", "Lambdat", "Zt", "RX", "X", 
        "RZX")), {
        CL <- solve(L, solve(L, Lambdat %*% Zt %*% sqrtW, system = "P"), 
            system = "L")
        CR <- solve(t(RX), t(X) %*% sqrtW - crossprod(RZX, CL))
        # if (fullHatMatrix) 
        #     crossprod(CL) + crossprod(CR)
        # else colSums(CR^2) + colSums(CL^2)

        # H = crossprod(CL) + crossprod(CR)
        # compute residual degrees of freedom as if using H, but use only CL and CR
        rdf_from_matrices(CL, CR)
    })
    rdf
}

# rdf.merMod(model)
# rdf.trace(model)




# system.time(replicate(100, rdf.trace(model)))
# system.time(replicate(100, rdf.merMod(model)))




#' Shrinkage metric for eBayes
#'
#' Shrinkage metric for eBayes quantifying the amount of shrinkage that is applied to shrink the maximum likelihood residual variance to the empirical Bayes posterior estimate
#'
#' @param sigmaSq maximum likelihood residual variance for every gene 
#' @param s2.post empirical Bayes posterior estimate of residual variance for every gene
#' 
#' @description Evaluates the coefficient from the linear regression of \code{s2.post ~ sigmaSq}. When there is no shrinkage, this value is 1.  Values less than 1 indicate the amount of shrinkage. 
#'
#' @import stats
#' @export
shrinkageMetric = function( sigmaSq, s2.post){
  fit = lm( s2.post ~ sigmaSq)
  as.numeric(coef(fit)[2])
}


