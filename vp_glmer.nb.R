
library(lme4)
require(MASS)
library(insight)

set.seed(101)
nrnd = 20
dd <- expand.grid(f1 = factor(1:3),
                  f2 = LETTERS[1:2], g=as.character(1:nrnd), rep=1:10,
          KEEP.OUT.ATTRS=FALSE)
dd$off = rnorm(nrow(dd))
mu <- (with(dd, as.integer(f1))) + (model.matrix(~g, dd) %*% rnorm(nrnd, 0, 1))

dd$y <- rnbinom(nrow(dd), mu = exp(mu*2), size = .5)
mean(dd$y)


system.time({
fit1 <- glmer.nb(y ~ offset(off) + f2 + (1|f1) + (1|g), data=dd, verbose=TRUE)
})
getME(fit1, "glmer.nb.theta")

# extract variance components
res = get_variance(fit1)



vc = with(res, c(fixed = var.fixed, random = res$var.intercept, residual = var.residual))

vc / sum(vc)

calcVarPart.glmerMod = function(fit){
	# Extracts variances from random effects
	# Computes total variance from fixed effects, 
	#     Then apportions across fixed effects
	# But Residuals is wrong because doesn't model distribution
	#
	# Note that get_variance(fit1) is very close, but doesn't apportion across fixed effects
	# Plus it is very slow for large datasets
	vc = variancePartition:::getVarianceComponents(fit)

	# extract distribution-specific variance
	vc$Residuals = get_variance_residual(fit)

	vc = unlist(vc)

	# create fractions
	varFrac = vc / sum(vc)

	# remove ".(Intercept)" string
	names(varFrac) = gsub("\\.\\(Intercept\\)", "", names(varFrac))

	varFrac
}

calcVarPart.glmerMod(fit1)

insight:::.compute_variances(fit, component = "all", name_fun = "get_variance", name_full = "random effect variances", model_component=NULL)





fit1 <- glmer.nb(y ~ offset(off) + f2 + (1|f1) + (1|g), data=dd, verbose=TRUE)

fit2 <- glmm.nb( y ~ offset(off) + f2, random = list(~ 1|f1, ~ 1|g), data=dd)



cor(residuals(fit2, type="pearson"), residuals(fit1, type="pearson"))	




# https://github.com/nyiuab/NBZIMM
# adapted from PQL, adds estimate of NB theta in each iteration
# PQL is problematic for count data with mean < 5
library(NBZIMM)
system.time({
fit2 = glmm.nb( y ~ offset(off) + f2, random = list(~1|g, ~1|f1), data=dd)
})
cor(residuals(fit2, type="pearson"), residuals(fit1, type="pearson"))	



library(NBZIMM)
system.time({
fitz = glmm.zinb( y ~ f1*f2, random = ~1|g, data=dd)
})
cor(residuals(fitz, type="pearson"), residuals(fit1, type="pearson"))	



library(glmmTMB)
system.time({
fit3 <- glmmTMB(y ~ f1*f2 + (1|g), data=dd, family = nbinom2)
})
cor(residuals(fit3, type="pearson"), residuals(fit1, type="pearson"))	


library(MASS)
system.time({
fit4 <- glmmPQL(y ~ f1*f2, random = ~ 1|g, data=dd, family=negative.binomial(2))
})
cor(residuals(fit4, type="pearson"), residuals(fit1, type="pearson"))	


library(rpql)
system.time({
XMM <- unname(model.matrix(fit1))
ZMM <- getME(fit1,"mmList");
names(ZMM) <- "cluster"
fit5 = rpql(y = getME(fit1, "y"), X = XMM, Z = list(cluster=getME(fit1, "Z")), id = list(cluster=dd$g), family = poisson(), lambda = 0, pen.type = "lasso")
})
cor(residuals(fit4, type="pearson"), residuals(fit1, type="pearson"))	


library(gamm4)
system.time({
fit6 <- gamm4(y ~ f1*f2, random = ~ 1|g, data=droplevels(dd[1:1000,]), family=negative.binomial(0.5))
})
cor(residuals(fit6$mer, type="pearson"), residuals(fit1, type="pearson")[1:1000])	


library(glmm)
system.time({
fit6 <- glmm(y ~ f1*f2, random = ~ 1|g, data=droplevels(dd[1:1000,]), family=negative.binomial(0.5), varcomps.names="g")
})
cor(residuals(fit6$mer, type="pearson"), residuals(fit1, type="pearson")[1:1000])	



system.time({
fit.test = nlme::lme( y ~ f1*f2, random = ~1|g, data=dd)
})


system.time({
fit.test2 = lme4::lmer( y ~ f1*f2 + (1|g), data=dd)
})


Variance components for GLMM
# https://easystats.github.io/insight/reference/get_variance.html
# https://royalsocietypublishing.org/doi/pdf/10.1098/rsif.2017.0213
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
# generating the design matrix is slow with get_variances():
# 	https://github.com/easystats/insight/blob/88ac255dfa3f559edaacb91e9e8a74e43046b40e/R/compute_variances.R
# 

library(insight)

get_variances(fit4)


0.193322+	0.288841+	0.410726	+0.107112


# Gene expression with counts
##############################

library('variancePartition')
library('edgeR')
library('BiocParallel')
data(varPartDEdata)

# filter genes by number of counts
isexpr = rowSums(cpm(countMatrix)>0.1) >= 5

# Standard usage of limma/voom
geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )


vp1 = lapply(1:1000, function(i){
	y = t(geneExpr$counts[i,,drop=FALSE])
	off = with(geneExpr$samples, log(lib.size * norm.factors))

	fit = glmer.nb(y ~ offset(off) + (1|Disease) + (1|Individual), metadata)

	calcVarPart.glmerMod(fit)
})
vp1 = do.call(rbind, vp1)



# fit = lmer(log(y) / log(off) ~ (1|Disease) + (1|Individual), metadata)

# calcVarPart(fit)


vobj = voomWithDreamWeights(geneExpr, ~ 1, metadata)


vp2 = fitExtractVarPartModel( vobj[1:1000,], ~ (1|Disease) + (1|Individual), metadata)

par(mfrow=c(1,3))
sapply(1:3, function(i){
 	plot(vp1[,i], vp2[,i], main=i)
 	abline(0, 1, col="red")
 	cor(vp1[,i], vp2[,i], method="spearman")
 })


# 1) Do VP on the count data directly and using glmer.nb
# 2) Use glmm.nb (nlme backend) to apply VST using covariates






























