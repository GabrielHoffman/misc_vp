
a = 270
b = 274
i = which((fit_mm$df.residual >=a) & (fit_mm$df.residual <=b))
n = nrow(info)
RSS = apply(fit_mm$residuals, 1, function(x) sum(x^2))

# plot(sigSq_m, fit_mm$sigma^2*fit_mm$df.residual)
# abline(0, 1, col="red")

# sigSq_m = fit_mm$sigma^2 * fit_mm$df.residual


plot(density(RSS[i] / var_true[i]), main="fit")
x = seq(0.1, 1000, length.out=10000)

# mixture of chisq works with 
dmixchisq = function(x, df){

	res = lapply(df, function(dff){
		dchisq(x, dff)
		})
	res = do.call(rbind, res)

	colSums(res) / length(df)
}

mixfit = dmixchisq(x, df=fit_mm$df.residual[i])

lines(x, mixfit, col="red")

# estimate df of chisq by KL divergence
KL = function(a, b, eps=1e-8){
	a[a<eps] = eps
	b[b<eps] = eps

	sum(a*(log(a) - log(b)))
}
h = function(df){
	KL(mixfit, dchisq(x, df)) #+ KL(dchisq(x, df), mixfit)
}

maxvalue = x[which.max(mixfit)]
opt = optimize(h, c(maxvalue*0.9, maxvalue*1.1))
opt
lines(x, dchisq(x, df=opt$minimum), col="green")



# gamma approx
h = function(v){
	KL(mixfit, dgamma(x, v[1], v[2])) #+ KL(dchisq(x, df), mixfit)
}

maxvalue = x[which.max(mixfit)]
opt = optim(c(a=maxvalue, b=1),h)
opt
lines(x, dgamma(x, opt$par[1], opt$par[2]), col="orange")





# empirical chisq approximation
f = function(df){
	sum(dchisq(y, df, log=TRUE))
}
opt = optimize(f, c(1, 1000), maximum=TRUE)
x = seq(0, 100, length.out=1000)
lines(x, dchisq(x, df=opt$maximum), col="blue")





# library(momentchi2)

# x = seq(0.01, 10000, length.out=100000)

# library(numDeriv)
# f = function(x){
# 	hbe(fit_mm$df.residual, x)
# }
# lines(x, grad(f, x))




df = seq(20, 400, length.out=10)
y = c(sapply(df, function(x){
	rchisq(1000, df=x)
	}))
plot(density(y, from=0))

df = c(sapply(df, function(x){
	rep(x, 1000)
}))

mixfit = dmixchisq(x, df=df)

lines(x, mixfit, col="red")


v = sapply(1:10000, function(i){
	z = rnorm(n)
	# (t(z) %*% (diag(1, n) - H) %*% (diag(1, n) - H) %*% z)[1]
	tcrossprod(t(z) %*% (diag(1, n) - H))[1]
	})

plot(density(v))




df = sum(diag((diag(1, n) - H )%*% (diag(1, n) - H)))

mixfit = dmixchisq(x, df=df)

lines(x, mixfit, col="red")


f = function(fit){

	suppressPackageStartupMessages(library(lme4))

	# number of samples
	n = nrow(fit@frame)

	H = as.matrix(hatvalues(fit, fullHatMatrix=TRUE))

	v = rgamma(1, 5, 1)

	z = rnorm(n, 0, sd=sqrt(v)) / sqrt(v)
	value = tcrossprod(t(z) %*% (diag(1, n) - H))[1]

	RSS = tcrossprod(t(fit@resp[['y']]) %*% (diag(1, n) - H))[1]

	data.frame(value, 
		RSS,
		rdf 	= n - 2*sum(diag(H)) + sum(H*H),
		redf 	= n - sum(H*H), 
		# rdf2 	= estimate_redf(H),
		v 		= v)
}
 # ~ Age + (1|Tissue) + (1|Individual) + (1|Batch)
res = fitVarPartModel( geneExpr,~ Age + (1|Individual)  + (1|Tissue)  + (1|Batch), info, fxn = f, REML=TRUE, BPPARAM=SerialParam() )

res = do.call(rbind, res)


plot(density(res$value, from=0))
mixfit = dmixchisq(x, df=res$rdf)
lines(x, mixfit, col="red")

lm(value ~ 0 + v, res)


plot(density(res$RSS / var_true, from=0))
mixfit = dmixchisq(x, df=res$rdf)
lines(x, mixfit, col="red")

lm(RSS ~ 0 + var_true, res)




plot((diag(1, n) - H), (diag(1, n) - H)%*%(diag(1, n) - H))

# S = (I-H)(I-H)
# z^T S z
# expectation of quadratic forms
S = (diag(1, n) - H)%*%(diag(1, n) - H)
mu = sum(diag( S ))

# variance of quadratic forms
# when the S is idempotent, SS = S so v = 2*mu
v = 2*sum(diag( S %*% S ))

# When S is not idempotent, this is an underdispersed chisq
2*mu
v


a = mu^2/v
b = mu /v

plot(x, dgamma(x, a, b))


estimate_redf = function(H){

	# S = (I-H)(I-H)
	# z^T S z
	# expectation of quadratic forms
	S = (diag(1, n) - H)%*%(diag(1, n) - H)
	mu = sum(diag( S ))

	# variance of quadratic forms
	# when the S is idempotent, SS = S so v = 2*mu
	v = 2*sum(diag( S %*% S ))

	# When S is not idempotent, this is an under/over-dispersed chisq
	# 2*mu
	# v

	# get moments of gamma		
	a = mu^2/v
	b = mu /v

	x = seq(1e-4, 10*mu, length.out=10000)

	# estimate df of chisq by KL divergence
	KL = function(a, b, eps=1e-8){
		a[a<eps] = eps
		b[b<eps] = eps

		sum(a*(log(a) - log(b)))
	}
	h = function(df){
		KL(y, dchisq(x, df)) + KL(dchisq(x, df), y)
	}
	y = dgamma(x, a, b)

	maxvalue = x[which.max(y)]
	opt = optimize(h, c(maxvalue*0.9, maxvalue*1.1))
	opt$minimum
}
estimate_redf(H)






# y = rgamma(x, a, b)
# mean(y)
# var(y)

























