#
# `b` and `d` are Gamma shape parameters and
# `a` and `c` are scale parameters.
# (All, therefore, must be positive.)
#
# Code from https://stats.stackexchange.com/questions/11646/kullback-leibler-divergence-between-two-gamma-distributions
KL.gamma <- function(a,b,c,d) {
  i <- function(a,b,c,d)
    - c * d / a - b * log(a) - lgamma(b) + (b-1)*(psigamma(d) + log(c))
  i(c,d,c,d) - i(a,b,c,d)
}

estimate_chisq_rdf = function(H){

	# estimate parameters for scaled chisq distribution
	a = sum(H*H) / tr(H)
	b = tr(H)^2 / sum(H*H)

	# estimate parameters for gamma 
	shape = b/2
	rate = 1/(2*a)
	scale = 1/rate

	# The KL divergence between a gamma distribution and a chisq distribution
	f = function(value){
		KL.gamma(value/2, 2, shape, scale) + KL.gamma(shape, scale, value/2, 2)
	}

	# Approximate the gamma distribution with a chisq 
	# estimate the degres of freedom for a chisq distribution that 
	# minimize the symmetric KL divergence
	fobj = optimize(f, c(1e-5, nrow(H)))
	df = fobj$minimum

	df
}





# eigen values of H
lambda = eigen(H)$values
lambda = lambda[lambda>1e-10]

# Monte Carlo sampling
obs_value = sapply( 1:100000, function(i){
	sum(rchisq(length(lambda), 1) * lambda)
})

est = mledist(obs_value, "gamma", fix.arg = list(rate = 1/2))

df = 2*est$estimate

plot(density(obs_value), main="mix chisq")

a = sum(H*H) / tr(H)
b = tr(H)^2 / sum(H*H)

# scaled chisq approx for Monte Carlo
values = a * rchisq(100000, b)
lines(density(values), col="blue")

# gamma PDF
x = seq(1e-4, max(obs_value), length.out=9999)
lines(x, dgamma(x, b/2, 1/(2*a)), col="green", lty=2)

# chisq approximation
lines(x, dchisq(x, tr(H)), col="red", lty=2)

df = estimate_chisq_rdf(H)
lines(x, dchisq(x, df), col="magenta")


n - sum(lambda)^2 / sum(lambda^2)


n - tr(H)^2 / sum(H*H)


satterthwaite_redf = function(H){
	n = nrow(H)
	S = (diag(1,n) - H) %*% (diag(1,n) - H)

	lambda = eigen(S)$values
	lambda = lambda[lambda>1e-10]

	sum(lambda)^2 / sum(lambda^2)
}

sum(lambda)
n - 2*tr(H) + sum(H*H)


n - .25*tr(H) + 0.5


 sum(lambda^2)

sum(eigen(S%*%S)$values)

A = diag(1,n) - 4*H + 6*H%*%H - 4*H%*%H%*%H + H%*%H%*%H%*%H 
tr(A)


n - 4*tr(H) + 6*sum(H*H) - 4*sum(H*(H%*%H)) + sum((H%*%H)*(H%*%H))



system.time(replicate(100, eigen(H, symmetric=TRUE, only.values=TRUE)$values))


system.time(replicate(100, n - 4*tr(H) + 6*sum(H*H) - 4*sum(H*(H%*%H)) + sum((H%*%H)*(H%*%H))))


system.time(replicate(100, partial_eigen(H, n=31)$values))

system.time(replicate(100, eigs_sym(H, k=31, opts = list(retvec = FALSE))$values))


system.time(replicate(100, satterthwaite_redf_fast(H, k=30)))


system.time(replicate(1000, rdf.satterthwaite(fit)))


eigen(H)$values[1:4]


eigen(H%*%H)$values[1:4]


ev = eigen(H)$values
# ev = ev[ev>1e-10]

rev(rep(1, n) - ev)

eigen(diag(1,n) - H)$values

sum(lambda)


S = (diag(1,n) - H) %*% (diag(1,n) - H)
lambda = eigen(S)$values
# lambda = lambda[lambda>1e-10]

satterthwaite_redf_fast = function(H, k){
	if( missing(k)) k = nrow(H)
	ev = eigs(H, k=k, opts = list(retvec = FALSE))$values

	lambda = rep(1, n)
	lambda[1:length(ev)] = 1 - ev
	lambda = lambda^2
	sum(lambda)^2 / sum(lambda^2)
}


H = as.matrix(hatvalues(fit, fullHatMatrix=TRUE))
		
satterthwaite_redf_fast(H)


satterthwaite_redf_fast(H)


ncol(getME(fit, "Z")) + ncol(getME(fit, "X"))


fit@Gp[length(fit@Gp)] + ncol(getME(fit, "X"))


#' Effective residual degrees of freedom
#'
#' Effective residual degrees of freedom for a linear mixed model using the Welch-Satterthwaite approximation
#' 
#' @param fit linear mixed model from \code{lmer()}
#' 
#' @return Effective residual degrees of freedom
#' 
#' @details
#' Compute the effective residual degrees of freedom, rdf, where the residual sum of squares (RSS) is approximately distributed as RSS/sigma^2 ~ \chi^2_{rdf}, where sigma^2 is the true residual variance.
#'
#' The residual degrees of freedom in a linear model is N - p - 1 where N is the sample size and p is the number of predictors.  In this case, rdf is defined based on the hat matrix H = X(X^TX)^{-1}X^T where rdf = n - tr(H)
#'
#'
#' @import lme4
#' @importFrom RSpectra eigs
#'
rdf.satterthwaite = function(fit){

	# number of non-zero eigen values of the hat matrix
	# is <= # random effect ncol + # fixed effect ncol
	k = fit@Gp[length(fit@Gp)] + ncol(lme4::getME(fit, "X"))

	# get hat matrix from linear mixed model
	H = hatvalues(fit, fullHatMatrix=TRUE)		

	# compute the first k eigen-values of H 
	# all subsequent eigen-values are zero 
	ev = RSpectra::eigs(H, k=k, opts = list(retvec = FALSE, tol=1e-4))$values
	
	# number of samples
	n = nrow(H)

	# compute lambda, the eigen-values of (I-H)(I-H)
	lambda = rep(1, n)
	lambda[1:length(ev)] = 1 - ev
	lambda = lambda^2

	# Welch-Satterthwaite approximation of degrees of freedom for 
	# chi-square
	sum(lambda)^2 / sum(lambda^2)
}
 # rdf.satterthwaite(fit)

		# linear time approximation
		# follows Hastie and Tibshirani. Generalized Additive Models. 1990
		# p54 for redf
   		# p305 for fast approximation


Compare emprical density if RSS/sigma^2
mixture of chi-square, 
and rdf.satterthwaite


Distribution of RSS is a mixture of chisq weighted by eigen values of (I-H)(I-H).  In linear regression, this simplifies substantially to RSS/sigma^2 ~ \chi^2_{n-p-1}.  This simplication depends on H being idempotent where eigen-values are either 1 or 0.  But the shrinkage in linear mixed models produces eigen-values between 0 and 1.

Naively computing eigen values of H is O(n^3).  Since only the first k eigen-values are non-zero they can be computed in O(n^2k).  But there is a linear time method too.  Since H = crossprod(A) + crossprod(B) where A and B are both low rank and sparse.  Compute lambda_a and lambda_b fast then then compute lambda.

library(onlinePCA)


# Initial data set
n <- 100        
d <- 50
X <- matrix(runif(n*d),n,d)
xbar <- colMeans(X)
pca0 <- eigen(cov(X))

# New observation
newx <- runif(d)

# Recursive PCA with secular equations
xbar <- updateMean(xbar, newx, n)
pca <- secularRpca(pca0$values, pca0$vectors, newx, n)#, center = xbar)
     
pca$values[1:4]

X_rbind = rbind(X, newx) #- xbar

eigen(cov(X_rbind))$values[1:4]







# my example with SVD

n <- 100        
d <- 50
X <- scale(matrix(runif(n*d),n,d))

dcmp.full = svd(X)
dcmp = eigen(X[,-1])

lambda = rep(0, n)
lambda[1:length(dcmp$d)] = dcmp$d^2

U = matrix(0, n, n)
U[,1:ncol(dcmp$u)] = dcmp$u

pca <- secularRpca(lambda, t(U), X[,1,drop=FALSE], n=500)


pca$values[1:4]
dcmp.full$d[1:4]^2
dcmp$d[1:4]^2








rdf.gamma.lambda = function(lambda, method = c("mean", "mode")){

	a = sum(lambda^2) / sum(lambda)
	b = sum(lambda)^2 /sum(lambda^2)

	# get mean and mode of gamma distribution
	# convert to mode of a chi-square distribution
	df = switch( method, 
			"mean" = a * b, 
			"mode" = a*(b-1) + 2)
	df
}





rdf.gamma = function(fit, method = c("mode", "mean")){

	method = match.arg(method)

	H = hatvalues(fit, fullHatMatrix=TRUE) 

	lambda = eigen(H)$values  


	nrow(fit@frame) - rdf.gamma.lambda( lambda, method)
}


rdf.merMod(fit, method="quad")
rdf.gamma(fit)





rdf.gamma = function(fit, method = c('1','2')){
	#
	# `b` and `d` are Gamma shape parameters and
	# `a` and `c` are scale parameters.
	# (All, therefore, must be positive.)
	#
	# Code from https://stats.stackexchange.com/questions/11646/kullback-leibler-divergence-between-two-gamma-distributions
	KL.gamma <- function(a,b,c,d) {
	  i <- function(a,b,c,d)
	    - c * d / a - b * log(a) - lgamma(b) + (b-1)*(psigamma(d) + log(c))
	  i(c,d,c,d) - i(a,b,c,d)
	}

	# Define trace function	
	tr = function(A) sum(diag(A))


	H = hatvalues(fit, fullHatMatrix=TRUE) 
	# estimate parameters for scaled chisq distribution
	a = sum(H*H) / tr(H)
	b = tr(H)^2 / sum(H*H)

	# lambda = eigen(H)$values  
	# lambda = lambda[lambda > 1e-10]
	# a = sum(lambda^2) / sum(lambda)
	# b = sum(lambda)^2 /sum(lambda^2)


	# estimate parameters for gamma 
	shape = b/2
	rate = 1/(2*a)
	scale = 1/rate

	# The KL divergence between a gamma distribution and a chisq distribution
	f = function(value){
		switch(method, 	'1' = KL.gamma(value/2, 2, shape, scale), 
						'2' = KL.gamma(shape, scale, value/2, 2))
	}

	# Approximate the gamma distribution with a chisq 
	# estimate the degres of freedom for a chisq distribution that 
	# minimize the symmetric KL divergence
	fobj = optimize(f, c(1e-5, 10*sum(lambda)))
	df = fobj$minimum

	df

}

method = "1"
x = seq(0, 60, length.out=1000)
y = sapply(x,f)
plot(x,log10(y))
x[which.min(y)]







