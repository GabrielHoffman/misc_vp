#  May 2020
# heritability from summary statstics

n = 100

beta = 0.01
x = scale(rnorm(n))
y = scale(x*beta + rnorm(n))

fit = lm(y~x)

se = coef(summary(fit))[2,2] * sqrt(n-1)


tstat = coef(summary(lm(y~x)))[2,3]

# is this true in expectation, rather than exact?
tstat^2 / (tstat^2 + n)

summary(fit)$r.squared


coef(summary(fit))




total heritiability should be

g = tstat^T Sigma^{-1} tstat

g/(g+n)


library(Rfast)
library(decorrelate)
library(aod)
library(Matrix)
n = 20003
p = 432

Sigma = bdiag(autocorr.mat(p/4, .9), autocorr.mat(p/4, .9), autocorr.mat(p/4, .9), autocorr.mat(p/4, .9))

beta = rnorm(p, 0, .01)
# beta[] = .02
X = scale(rmvnorm(n, rep(0, p), drop0(Sigma, 0.01)))
eta = X%*%beta
y = eta + rnorm(n)
y_scale = scale(y)

var(eta) / var(y)

df = lapply(1:ncol(X), function(i){
	fit = lm(y~X[,i])
	tstat = coef(summary(fit))[2,3]
	data.frame(rsq = summary(fit)$r.squared,
		r.sq.t = tstat^2 / (tstat^2 + n), beta=coef(fit)[2])
})
df = do.call(rbind, df)



# Pasaniuc: univariate beta
a = t(df$beta) %*% solve(cora(X) + diag(0, p), df$beta)
(n*a - p) / (n-p)

ecl = eclairs(X, compute="cor")
k = tr_cp(ecl)
# getCov(ecl) %*% beta
(n*quadForm(ecl, df$beta) -k) / (n-k)

# W = crossprod(decorrelate(decorrelate(df$beta, ecl), ecl))*n
W = crossprod(decorrelate(df$beta, ecl))*n

W/(W+n)




plot(beta, decorrelate(decorrelate(df$beta, ecl), ecl))





plot(solve(chol(cora(X)), df$beta))



plot(cora(X) %*% beta, df$beta)





plot( beta, solve(cora(X), df$beta))






plot(beta, solve(chol(cora(X)), df$beta))


plot(df$beta, decorrelate(coef(fit)[-1], ecl))





fit = lm(y~X)
summary(fit)$r.squared



X = scale(rmvnorm(n, rep(0, p), drop0(Sigma, 0.01)))
sum(diag(solve(crossprod(X), crossprod(X))))

ecl = eclairs(X)
sum(diag(solve(getCor(ecl), crossprod(X)))) / n

tr_cp(ecl)

# trace(solve(eclairs(X), crossprod(X)))
tr_cp = function(ecl){

	with(ecl, sum(dSq/((1-lambda)*dSq + lambda*nu)))
}


# tstat = coef(summary(fit))[2:(p+1),3]

# g = tstat^T Sigma^{-1} tstat

# g = crossprod(tstat, solve(Sigma, tstat))
# g = crossprod(tstat, Sigma %*% tstat)
# g = crossprod(tstat)
# g/(g+n)


# plug wald statistic into heritability 
res = wald.test(vcov(fit), coef(fit), 1:(p+1))
g = res$result$chi2[1]
g/(g+n)


g = t(coef(fit)) %*% solve(vcov(fit)) %*% coef(fit)
g/(g+n)





sSq = sum(fit$residuals^2) / (n-p)

# tstat = coef(summary(fit))[,3]
# beta = coef(fit)

# g = t(beta) %*% crossprod(cbind(1, X)) %*% beta / sSq
g = crossprod(X%*% beta) / sSq
g/(g+n)


Sigma = solve(crossprod(cbind(1, X))) * sSq
g = t(coef(fit)) %*% solve(Sigma) %*% coef(fit)
g/(g+n)


sum(solve(crossprod(X)/(n-1)) %*%  beta[-1])
ecl = eclairs(X, lambda=0)
getCov(ecl) %*%  beta[-1]

quadForm(ecl, beta[-1])











