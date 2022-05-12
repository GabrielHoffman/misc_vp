

library(Bolstad)
library(Rfast)
library(parallel)
library(mvtnorm)

library(data.table)
library(ggplot2)
library(ashr)

p = 50
n = 100

covariate = rnorm(n)

beta = c(1, 1, 1, rep(0, p-3))

Y = lapply(1:p, function(i){
	covariate * beta[i] + rnorm(n)*.5
	})
Y = do.call(rbind, Y)


evalfit = function(params, extra=FALSE){
	# b0 = params[1:2]
	# V0 = matrix(params[3:6], ncol=2)
	b0 = c(0,0)
	V0 = diag(c(params[1], params[2]))
	fitList = lapply( 1:p, function(j){
		info = data.frame(y = Y[j,], x = covariate)

		bayes.lm( y ~ x, info, prior = list(b0 = b0, V0 = V0))
	})
	
	ll_data = sum(sapply(fitList, logLik)) 
	ll_prior = sum(sapply(fitList, function(fit) dmvnorm( coef(fit), b0, V0, log=TRUE)))

	# res = ll_data + ll_prior
	res = sum(dmvnorm( beta.mle, b0, V0, log=TRUE))

	if(extra){
		beta = do.call(rbind, lapply(fitList, coef))
		se = do.call(rbind, lapply(fitList, function(fit) 
			sqrt(diag(fit$post.var))))

		negProb = pnorm(0, beta[,2], se[,2])

		lfsr = sapply(negProb, function(value) ashr::compute_lfsr( value, 0))
		sum(lfsr < 0.05)


		attr(res, "beta") = beta
		attr(res, "se") = se
		attr(res, "lfsr") = lfsr
	}
	res
}


beta.mle = lapply( 1:p, function(j){

	info = data.frame(y = Y[j,], x = covariate)

	coef(lm( y ~ x, info))
})
beta.mle = do.call(rbind, beta.mle)


opt1 = optim( c(1, 1), evalfit, 
	lower = 1e-8, upper=1e6, control=list(fnscale=-1), method='L-BFGS-B' )

opt2 = list(par = apply(beta.mle, 2, var))
res1 = evalfit( opt1$par, extra=TRUE)
res2 = evalfit( opt2$par, extra=TRUE)

beta1 = attr(res1, "beta")
colnames(beta1) = colnames(beta.mle)
beta2 = attr(res2, "beta")
colnames(beta2) = colnames(beta.mle)


df = rbind(rbind( data.frame(Method = 'beta.mle', beta.mle, beta.true=beta),
	data.frame(Method = 'EB1', beta1, beta.true=beta)),
	data.frame(Method = 'EB2', beta2, beta.true=beta) )
df = data.table(df)

# ggplot(df, aes(x, color=Method)) + geom_density() + theme_classic()

df[,sum((x-beta.true)^2), by='Method']

ggplot(df, aes(beta.true, x, color=Method))  + theme_classic() + geom_point() + facet_grid(~Method)






j=1
b0 = c(0,0)
V0 = diag(c(1e6, 1e6))
info = data.frame(y = Y[j,], x = covariate)
fit = bayes.lm( y ~ x, info, prior = list(b0 = b0, V0 = V0))
vcov(fit)
summary(lm( y ~ x, info))


value = seq(-8, 5, length.out=100)
score = mclapply( value, function(x){
	V0 = diag(c(10^x, 10^x))
	fit = bayes.lm( y ~ x, info, prior = list(b0 = b0, V0 = V0))
	logLik(fit) 
	dmvnorm( coef(fit), b0, V0, log=TRUE)
	}, mc.cores=5)
score = unlist(score)

x = value[which.max(score)]

plot(value, score)




hist(attr(res1, "lfsr"))
sum(attr(res1, "lfsr") < 0.05)





value = seq(-4, 5, length.out=100)
score = mclapply( value, function(x){
	evalfit( c(10^x, 10^x))
	}, mc.cores=5)
score = unlist(score)

x = value[which.max(score)]

plot(value, score)







fit = lm(t(Y) ~ covariate)
tstat = sapply( summary(fit), function(x) coef(x)[2,'t value'] )

res.ash = ash.workhorse(tstat, rep(1, p), alpha=1, pointmass=FALSE)

sum(res.ash$result$lfsr < 0.05)








evalfit = function(params, extra=FALSE){
	b0 = c(0)
	V0 = (params
	fitList = lapply( 1:p, function(j){
		info = data.frame(tstat = tstat, covariate)

		bayes.lm( tstat ~ covariate, info, prior = list(b0 = b0, V0 = V0))
	})
	
	res = sum(sapply(fitList, logLik))
	if(extra){
		beta = do.call(rbind, lapply(fitList, coef))
		se = do.call(rbind, lapply(fitList, function(fit) sqrt(diag(fit$post.var))))

		negProb = pnorm(0, beta[,2], se[,2])

		lfsr = sapply(negProb, function(value) ashr::compute_lfsr( value, 0))

		attr(res, "beta") = beta
		attr(res, "se") = se
		attr(res, "lfsr") = lfsr
	}
	res
}


p = 10000
x = rnorm(p)
y = x + rnorm(p) + 10

coef(lm(y ~ x))[2]
coef(lm(x ~ y))[2]
mean(x-y)

coef(lm(y ~ 0+x))[1]
coef(lm(x ~ 0+y))[1]

m = 96
reg = 11
M = matrix(runif(96*11), m, 11)

sum(apply(M, 1, min) < 0.05)





library(lavaan)

n = 50
x = sample(0:2, n, replace=TRUE)
y1 = x + rnorm(n)
y2 = y1 + rnorm(n)


df = data.frame(SNP = x, Y1obs = y1, Y2obs = y2)
cor(df)



model = " 	
			# latent variables
			y1 =~ SNP
			y2 =~ y1

			# regression
			Y1obs ~ y1
			Y2obs ~ y2
		"

fit = sem(model, df)
BIC(fit)



df2 = data.frame(SNP = x, Y1obs = y2, Y2obs = y1)

fit2 = sem(model, df2)
BIC(fit2)



fitSEM = function(SNP, Y1, Y2){

	Data <- data.frame(X = SNP, Y = Y1, M = Y2)
	model <- ' # direct effect
	             Y ~ c*X
	           # mediator
	             M ~ a*X
	             Y ~ b*M
	           # indirect effect (a*b)
	             ab := a*b
	           # total effect
	             total := c + (a*b)
	         '
	fit <- sem(model, data = Data)
}



fit1 = fitSEM( x, y1, y2)
fit2 = fitSEM( x, y2, y1)

BIC(fit1)
BIC(fit2)







