
f = function(fit){

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

		x = seq(1e-4, 10*mu, length.out=10000)

		# get moments of gamma		
		a = mu^2/v
		b = mu /v

		# pdf of gamma
		y = dgamma(x, a, b)

		# estimate df of chisq by KL divergence
		KL = function(a, b, eps=1e-8){
			a[a<eps] = eps
			b[b<eps] = eps

			sum(a*(log(a) - log(b)))
		}
		h = function(df){
			KL(y, dchisq(x, df)) + KL(dchisq(x, df), y)
		}

		maxvalue = x[which.max(y)]
		opt = optimize(h, c(maxvalue*0.5, maxvalue*2))
		data.frame(rdf = opt$minimum, a=a, b=b)
	}

	# Generate samples from mixture of chisq
	# 	this requires eigen decomp
	# but can use single scaled chisq for speed
	# DOI:10.1348/000711009X449771
	estimate_rdf_mc = function(H){

		# eigen values of H
		lambda = eigen(H)$values
		lambda = lambda[lambda>1e-10]

		# Monte Carlo sampling
		obs_value = sapply( 1:100000, function(i){
			sum(rchisq(length(lambda), 1) * lambda)
		})

		est = mledist(obs_value, "gamma", fix.arg = list(rate = 1/2))

		df = 2*est$estimate
		df
	}

	suppressPackageStartupMessages(library(lme4))
	suppressPackageStartupMessages(library(fitdistrplus))

	# number of samples
	n = nrow(fit@frame)

	H = as.matrix(hatvalues(fit, fullHatMatrix=TRUE))

	# sum(residuals(fit)^2)
	RSS = tcrossprod(t(fit@resp[['y']]) %*% (diag(1, n) - H))[1]

	r = fit@resp[['y']] - fitted(fit)
	sum(r^2)

	res = estimate_redf(H)

	# fixed effect model
	fit_fixed = lm(fit@frame[,1] ~ ., fit@frame[,-1])

	data.frame(	RSS 	= RSS, 
				RSS_fixed = sum(residuals(fit_fixed)^2),
				rdf 	= n - 2*sum(diag(H)) + sum(H*H), 
				rdf2 	= res$rdf, 
				rdf.mc 	= n - estimate_rdf_mc(H),
				rdf.fixed = n - sum(hatvalues(fit_fixed)),
				a 		= res$a, 
				b 		= res$b)
}

 # ~ Age + (1|Tissue) + (1|Individual) + (1|Batch)
res = fitVarPartModel( geneExpr, ~ Age + (1|Individual)  + (1|Tissue) + (1|Batch), info, fxn = f, REML=FALSE, BPPARAM=SnowParam(6) )

res = do.call(rbind, res)

# fixed effect model
plot(density(res$RSS_fixed / var_true), main="fixed", lwd=2)
x = seq(0, 1000, length.out=10000)
lines(x, dchisq(x, df=res$rdf.fixed), col="red")


# mixed model
plot(density(res$RSS / var_true), main="mixed", lwd=2, ylim=c(0, 0.03))
x = seq(0, 1000, length.out=1000)
lines(x, dchisq(x, df=mean(res$rdf2)), col="red")


# mixture of chisq works with 
dmixchisq = function(x, df){

	res = lapply(df, function(dff){
		dchisq(x, dff)
		})
	res = do.call(rbind, res)

	colSums(res) / length(df)
}


value = dmixchisq(x, res$rdf)
lines(x, value, col="green")

value = dmixchisq(x, res$rdf.mc)
lines(x, value, col="blue")


# mixture of gamma
dmixgamma = function(x, a,b){

	res = lapply(1:length(a), function(i){
		dgamma(x, a[i],b[i])
		})
	res = do.call(rbind, res)

	colSums(res) / length(a)
}

value = dmixgamma(x, res$a, res$b)
lines(x, value, col="blue")


coef(lm(RSS_fixed ~ 0 + var_true, res))
coef(lm(RSS ~ 0 + var_true, res))



# Approximate a mixture of chisq's with a single chisq
# following DOI:10.1348/000711009X449771
# H = diag(1, 10)
# diag(H) = c(rep(1, 5), .6, .7, rep(0, 9))
lambda = eigen(H)$values
lambda = lambda[lambda>1e-10]

tr(H)
tr(H%*%H)

# lambda = c(1,1,1) 

tr = function(A) sum(diag(A))

c = tr(H) / length(lambda)

obs_value = sapply( 1:100000, function(i){
	sum(rchisq(length(lambda), 1) * lambda)
	})

plot(density(obs_value, from=0))

# gamma estimate
est = mledist(obs_value, "gamma")
x = seq(1e-4, max(obs_value), length.out=9999)
lines(x, dgamma(x, est$estimate[1], est$estimate[2]), col="red", lty=2)
est$estimate

# chisq estimate
est = mledist(obs_value, "gamma", fix.arg = list(rate = 1/2))
x = seq(1e-4, max(obs_value), length.out=9999)
lines(x, dgamma(x, est$estimate[1], 1/2), col="blue", lty=2)

# df
2*est$estimate
tr(H)
tr(H%*%H)

mean(obs_value)
var(obs_value)


# Generate samples from mixture of chisq
# 	this requires eigen decomp
# but can use single scaled chisq for speed
# DOI:10.1348/000711009X449771
estimate_rdf_mc = function(H){

	# eigen values of H
	lambda = eigen(H)$values
	lambda = lambda[lambda>1e-10]

	# Monte Carlo sampling
	obs_value = sapply( 1:100000, function(i){
		sum(rchisq(length(lambda), 1) * lambda)
	})

	est = mledist(obs_value, "gamma", fix.arg = list(rate = 1/2))

	df = 2*est$estimate
	df
}

n - estimate_rdf_mc(H)

n - 2*sum(diag(H)) + sum(H*H)



x = seq(1e-4, max(obs_value), length.out=9999)

df = sum(lambda)
lines(x, dchisq(x, df), col="red", lty=2)

# approx 1
df = length(lambda)
lines(x, dchisq(x/c, df), col="blue", lty=2)


# approx 2
a = sum(lambda^2) / sum(lambda)
b = sum(lambda)^2 / sum(lambda^2)
lines(x, dchisq(x/a, b, ncp=a), col="green", lty=2)


c = 7
obs_value = c* rchisq(10000, 1)

mean(obs_value)
var(obs_value)

# this gives mean of c and var of 2*c^2
plot(density(obs_value, from=0))

x = seq(1e-4, max(obs_value), length.out=9999)
lines(x, dchisq(x, c, ncp=c), col="red", lty=2)



# mixture of chisq works with 
dmixchisq = function(x, lambda){

	res = lapply(lambda, function(lamb){
		lamb* dchisq(x, 1)
		})
	res = do.call(rbind, res)

	colSums(res) / length(lambda)
}


mixfit = dmixchisq(x, lambda)

lines(x, mixfit, col="red")













