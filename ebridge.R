# Gabriel Hoffman
# May 11, 2022

#' Empirical Bayes ridge regression
#'
#' Empirical Bayes ridge regression with shrinks covariance
#'
#' @param y vector of matrix of response variables
#' @param X design matrix
#' @param lambda shrinkage parameter.  If \code{NULL}, estimated by empirical Bayes.
#' @param k rank of \code{eclairs()} decomposition
#'
#' @details
#' OLS:  \eqn{\beta = (X^TX)^{-1} X^T y}
#' ridge:  \eqn{\beta = (X^TX + I\lambda)^{-1} X^T y}
#' EB ridge: \eqn{\beta = (X^TX*(1-\lambda) + I\lambda)^{-1} X^T y}
ebridge = function(y, X, k=min(dim(X)), lambda=NULL){

	# get whitening transformation
	ecl = eclairs(X, lambda=lambda, compute="corr")

	# whiten X, and transform y
	beta = crossprod(decorrelate(X, ecl, alpha=-1), y) 

	# how do I get covariance, and standard errors?

	# rescale coefs 
	# since decorrelate() standardizes the variances of X
	list( coefficients = beta / nrow(X), 
		lambda = ecl$lambda)
}




# EB ridge regression
library(Rfast)
library(decorrelate)
library(MASS)
library(Matrix)

n = 1003
p = 535

Sigma = drop0(autocorr.mat(p, .9), 0.01)

X = rmvnorm(n,rep(0,p), Sigma)
beta.true = rnorm(p)
y = X %*% beta.true + rnorm(n)

fit = lm.ridge(y ~ 0 + ., data.frame(X), lambda=10^seq(-5, 6, length.out=199))
# plot(fit$GCV)

 fit$lambda[which.min(fit$GCV)]

beta.ridge = fit$coef[,which.min(fit$GCV)]

fit = ebridge(y, X)
fit$lambda

beta.eb = c(fit$coefficients)

par(mfrow=c(1:2))
plot(beta.true, beta.ridge)
abline(0,1, col="red")
plot(beta.true, beta.eb)
abline(0,1, col="red")


cor(beta.true, beta.ridge)
cor(beta.true, beta.eb)









