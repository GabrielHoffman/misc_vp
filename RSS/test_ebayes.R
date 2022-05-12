

# FOR LMM 
# H and I-H are not idempotent

library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(lmerTest)
library(fitdistrplus)

data(varPartData)

# info = rbind(info, info, info)
lvl = levels(info$Individual)
info = droplevels(info[info$Individual %in% lvl[1:4],])
n_genes = 1000
var_true = rgamma(n_genes, 5, 1)
# var_true[] = 10

geneExpr = lapply(1:n_genes, function(j){

	beta_ID = rnorm(nlevels(info$Individual))
	beta_Tissue = rnorm(nlevels(info$Tissue))
	beta_Batch = rnorm(nlevels(info$Batch))

	y = model.matrix(~0+Individual, info) %*% beta_ID + 	
		# model.matrix(~0+Tissue, info) %*% beta_Tissue + 
		# model.matrix(~0+Batch, info) %*% beta_Batch +
		rnorm(nrow(info), 0, sd=sqrt(var_true[j])) 
	t(y)
	})
geneExpr = do.call(rbind, geneExpr)

form = ~ Age + (1|Individual) 
# form = ~ Age + (1|Individual)  + (1|Tissue)  + (1|Batch)

# fit = dream( geneExpr, ~ Age + Individual+ Batch + Tissue, info)
# dsgn = model.matrix(~ Age + Individual + Batch + Tissue, info)
dsgn = model.matrix(subbars(form), info)
fit = lmFit( geneExpr, dsgn)

# linear mixed model
fit_mm = dream( geneExpr, form, info, BPPARAM=SnowParam(6))

i = 129

y = t(t(geneExpr[i,]))
form2 = as.formula(paste0('y ~ ', as.character(form)[-1], sep='', collapse=''))
it = lmer(form2, info, REML=TRUE)

H = as.matrix(hatvalues(it, fullHatMatrix=TRUE))
		
plot( residuals(it), crossprod(y, diag(1,nrow(H))-H))
plot( residuals(fit_mm)[129,], crossprod(y, diag(1,nrow(H))-H))

# it2 = lm(y ~ Age + Tissue + Individual + Batch, info)

# summary(it)

# summary(it2)

# H = hatvalues(it, fullHatMatrix=TRUE)

# sum(diag(crossprod(diag(1, nrow(H)) - H)))
# n = nrow(H)
# n - 2*sum(diag(H)) + sum(H*H)

# n - 1.25*sum(diag(H)) + 0.5


# system.time(replicate(1000, sum(diag(crossprod(H)))))
# system.time(replicate(1000, sum(H*H)))


# system.time(replicate(1000, n - 2*sum(diag(H)) + sum(H*H)))

# system.time(replicate(1000, n - 1.25*sum(diag(H)) + 0.5))


f = function(fit){

	suppressPackageStartupMessages(library(lme4))
	tr = function(A) sum(diag(A))

	satterthwaite_redf = function(H){
		n = nrow(H)
		S = (diag(1,n) - H) %*% (diag(1,n) - H)

		lambda = eigen(S)$values
		lambda = lambda[lambda>1e-10]

		sum(lambda)^2 / sum(lambda^2)
	}

	# number of samples
	n = nrow(fit@frame)

	if( n < 200){
		# more accurate but quadratic time
		H = as.matrix(hatvalues(fit, fullHatMatrix=TRUE))
		# n - 2*sum(diag(H)) + sum(diag(crossprod(H)))
		# redf = n - 2*sum(diag(H)) + sum(H*H)

		# Satterwait
		redf = satterthwaite_redf(H)
	}else{
		# linear time approximation
		# follows Hastie and Tibshirani. Generalized Additive Models. 1990
		# p54 for redf
   		# p305 for fast approximation
		h.diag = hatvalues(fit)
		redf = n - 1.25*sum(h.diag) + 0.5
	}
	redf
}

f = rdf.satterthwaite


 # ~ Age + (1|Tissue) + (1|Individual) + (1|Batch)
res = fitVarPartModel( geneExpr, form, info, fxn = f, REML=TRUE, BPPARAM=SnowParam(6) )

df.values = unlist(res)

df.residual.orig = fit_mm$df.residual 

fit_mm$df.residual = df.values

#### check discrepency between residuals and sigmasq


# fit$sigma^2 is already corrected by n-1
n = nrow(info)
rdf = n - ncol(fit$design)
RSS = apply(residuals(fit, geneExpr), 1, function(x) sum(x^2))
# plot(sigSq, fit$sigma^2 * df); abline(0, 1, col="red")
# sigSq = fit$sigma^2 * df
plot(density(RSS / var_true, from=0), main="fit")
x = seq(0, 1000, length.out=10000)
lines(x, dchisq(x, df=rdf), col="red")


a = 0
b = 1000
i = which((fit_mm$df.residual >=a) & (fit_mm$df.residual <=b))
n = nrow(info)
RSS = apply(fit_mm$residuals, 1, function(x) sum(x^2))
plot(density( RSS[i] / var_true[i], from=0), main="fit")
x = seq(0, 1000, length.out=10000)


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

y = RSS / var_true
# empirical chisq approximation
f = function(df){
	sum(dchisq(y[i], df, log=TRUE))
}
opt = optimize(f, c(1, 10000), maximum=TRUE)
x = seq(0, 1000, length.out=10000)
lines(x, dchisq(x, df=opt$maximum), col="blue")





# is RSS computed correctly? for lmm since 
# RSS = e^T (I - H)^T (I - H) e where e is the true error

# e^T (I - H) e = e^T e  - e^T H e

# use a mixture of chisquared distrubiotn

fit_mm_mod = fit_mm

RSS = apply(fit_mm$residuals, 1, function(x) sum(x^2))
sigSqm = RSS / fit_mm$df.residual 
fit_mm_mod$sigma = sqrt(sigSqm)

fit_eb = limma::eBayes(fit)
fit_eb_mm = limma::eBayes(fit_mm)
fit_eb_mm_mod = limma::eBayes(fit_mm_mod)

# Empirical Bayes results
par(mfrow=c(1,2))
plot(density((fit_eb$df.prior*fit_eb$s2.prior)/fit$sigma^2), main="fit")
x = seq(0, 100, length.out=1000)
lines(x, dchisq(x, fit_eb$df.prior), col="red")

plot(density((fit_eb_mm$df.prior*fit_eb_mm$s2.prior)/fit_eb_mm$sigma^2), main = "lmm")
x = seq(0, 100, length.out=1000)
lines(x, dchisq(x, fit_eb_mm$df.prior), col="red")

fit_eb$df.prior
fit_eb_mm$df.prior
fit_eb_mm_mod$df.prior

fit_eb$s2.prior
fit_eb_mm$s2.prior
fit_eb_mm_mod$s2.prior


plot(fit$sigma^2, fit_mm$sigma^2)
abline(0, 1, col='red')


plot(fit_eb$s2.post, fit_eb_mm$s2.post)
abline(0, 1, col='red')


par(mfrow=c(4,3))
hist(fit$sigma^2)
hist(fit$df.residual)
abline(v=fit$df.residual[1], col="red")
hist(fit_eb$s2.post)

hist(fit_mm$sigma^2)
hist(fit_mm$df.residual)
hist(fit_eb_mm$s2.post)

plot(fit_mm$sigma^2, fit_eb_mm$s2.post)
abline(0, 1, col="red")


d1 = density(fit$sigma^2, from=0)
d2 = density(fit_eb$s2.post,from=0)
ymax = max(c(max(d1$y), max(d2$y)))
plot(d1, col="blue", lty=2, ylim=c(0, ymax), main="Fixed effects")
lines(d2, col="green")

d1 = density(fit_mm$sigma^2, from=0)
d2 = density(fit_eb_mm$s2.post,from=0)
ymax = max(c(max(d1$y), max(d2$y)))
plot(d1, col="blue", lty=2, ylim=c(0, ymax), main="Mixed effecs")
lines(d2, col="green")

df_mse = data.frame(sigSq = mean((fit$sigma^2 - var_true)^2),
					sigSq_mm = mean((fit_mm$sigma^2 - var_true)^2),
					sigSq_post = mean((fit_eb$s2.post - var_true)^2),
					sigSq_mm_post = mean((fit_eb_mm$s2.post - var_true)^2), 
					sigSq_mm_mod_post = mean((fit_eb_mm_mod$s2.post - var_true)^2))

barplot(as.matrix(df_mse))

# CONCLUSION
s2.post for LMM reduces MSE

need to translate the new df.prior to topTable
currently, fit_eb_mm$p.value is already computed by dream
Need to adjust t-stat and then df.total = df. + df.prior

in limma the df is based on the residual df


in dream df is the Satterthwaite method
	how to combine df.prior with df.Satterthwaite

# Residual_effective_degrees_of_freedom
https://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics)#Residual_effective_degrees_of_freedom

Generalized Additive Models, Hastie and Tibshirani, 1990
  p54 for redf
   p305 for fast approximation
https://books.google.com/books?id=qa29r1Ze1coC&lpg=PR3&dq=Hastie%2C%20T.%20J.%2C%20and%20Tibshirani%2C%20R.%20J.%20(1990)%2C%20Generalized%20Additive%20Models%2C%20London%3A%20Chapman%20and%20Hall.&pg=PA54#v=onepage&q&f=false


fit_eb$df.total
fit_eb_mm$df.total

coef = "(Intercept)"
p = topTable(fit_eb, coef, number=Inf, sort.by="none")[,'P.Value']
sum(p < 0.05) / length(p)


p = topTable(fit_mm, coef, number=Inf, sort.by="none")[,'P.Value']
sum(p < 0.05) / length(p)

p = topTable(fit_eb_mm, coef, number=Inf, sort.by="none")[,'P.Value']
sum(p < 0.05) / length(p)

p.new = 2*pt(abs(fit_eb_mm$t[,coef]), df=df.residual.orig[,coef] + fit_eb_mm$df.prior, lower.tail=FALSE)
sum(p.new < 0.05) / length(p.new)


p.new = 2*pt(abs(fit_eb_mm$t[,coef]), df=df.residual.orig[,coef], lower.tail=FALSE)
sum(p.new < 0.05) / length(p.new)


p.new = 2*pt(abs(fit_eb_mm_mod$t[,coef]), df=df.residual.orig[,coef] + fit_eb_mm_mod$df.prior, lower.tail=FALSE)
sum(p.new < 0.05) / length(p.new)






