# Gabriel Hoffman
# Feb 10, 2021
# 
# Evalaute evidence for direction of causality

library(cit)

# Evalute forward and reverse model, alowing for covariates
eval_regulatory_direction = function(SNP, Y_OCR, Y_gene, covariates=NULL){

	# SNP -> Y_OCR -> Y_gene (forward)
	res1 = cit.cp( SNP, Y_OCR, Y_gene, covariates)

	# SNP -> Y_gene -> Y_OCR (reverse)
	res2 = cit.cp( SNP, Y_gene, Y_OCR, covariates)

	data.frame( p.forward = res1['p_cit'],
				p.reverse = res2['p_cit'])
}


# perform analysis on 100 simulated triplets 
df_res = lapply(1:100, function(i){

	n = 500
	alpha = rnorm(1)
	beta = rnorm(1)

	# Simulate for forward model
	SNP = sample(0:2, n, replace=TRUE)
	y_ocr = SNP*alpha + rnorm(n)
	y_gene = y_ocr*beta + rnorm(n)

	# estimate direction based on data
	eval_regulatory_direction( SNP, y_ocr, y_gene )
	})
df_res = do.call(rbind, df_res)

# Bonferroni cutoff
cutoff = 0.05 / nrow(df_res)

# Assign direction based on results and Bonferroni cutoff
df_res$direction = "Unknown"
df_res$direction[with(df_res, (p.forward < cutoff) & (p.reverse > cutoff))] = "Forward"
df_res$direction[with(df_res, (p.forward > cutoff) & (p.reverse < cutoff))] = "Reverse"
df_res$direction[with(df_res, (p.forward > cutoff) & (p.reverse > cutoff))] = "Independent"

