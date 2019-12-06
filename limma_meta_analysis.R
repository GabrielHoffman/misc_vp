# Gabriel Hoffman
# December 5, 2019

library(variancePartition)
library(metafor)

# perform meta-analysis given fits from limma
run_meta_analysis = function( fitList, coef, method="REML" ){

	if( length(coef) > 1){
		stop("Only supports one coefficient term")
	}

	# extract logFC and se for each model fit
	resList = lapply( fitList, function(fit){
		res = topTable(fit, coef, sort.by="none", number=Inf)
		res$se = res$logFC /  res$t
		res
	})

	# extract features names
	featureNames = lapply(resList, function(res){
		rownames(res)
		})
	featureNames = unique(unlist(featureNames))

	# check that each fit has all features
	testIdentical = lapply(resList, function(res){
		identical(sort(rownames(res)), sort(featureNames))
		})
	if(  any(!unlist(testIdentical)) ){
		stop("each entry in fitList must have all features with the same names")
	}

	# apply meta analysis for each feature
	resRMA = lapply(featureNames, function(feature){
		# extract summary statistics for this feature from all studies
		res = lapply(resList, function(res){
			res[feature,]
			})
		df = do.call("rbind", res)
		rma(yi=logFC, sei=se, data = df, method=method)
		})
	names(resRMA) = featureNames

	# extract results
	# Currently exatract simple results
	# other results may be relevant if meta-analysing more than two studies
	resTable = lapply(resRMA, function(x){
		data.frame(	beta 	= x$beta,
					se 		= x$se,
					pvalue 	= x$pval)
		})
	do.call("rbind", resTable)
}


# get data from example analysis for dream
data(varPartData)

# simulate a second dataset
geneExpr2 = matrix(rnorm(length(geneExpr)), nrow(geneExpr), ncol(geneExpr))
rownames(geneExpr2) = rownames(geneExpr)
colnames(geneExpr2) = colnames(geneExpr)

# Specify formula
form <- ~ Batch 

# run dream, here same as lmFit()
fit1 = dream( geneExpr, form, info)
fit2 = dream( geneExpr2, form, info)

# run eBayes
fit1 = eBayes(fit1)
fit2 = eBayes(fit2)


# create a list of dream or eBayes results to meta analyze
fitList = list(fit1, fit2)

# specifiy coefficient to test
coef = "Batch2"

# random effect meta-analysis
res_random = run_meta_analysis( fitList, coef, method="REML")

# fixed effect meta-analysis
res_fixed = run_meta_analysis( fitList, coef, method="FE")









