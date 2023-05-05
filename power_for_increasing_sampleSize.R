# Gabriel Hoffman
# May 3, 2023
#
# Evaluate power using Cohen's d and observed beta and se_beta


# Cohen's D
############
# Consider comparing two groups, A and B.  The power is related to 
# the difference in means between the groups (i.e. effect size) and the pooled standard deviation.
# Cohens d is the effect size divided by this pooled standard deviation


library(ggplot2)
library(pwr)

# correct for 15K genes for multiple testing
n_tests = 17128
n_subjects = c(44, 100, 200)

grid = expand.grid(delta = seq(0, 1, length.out=100), n_subjects = n_subjects)

res = lapply(seq(nrow(grid)), function(i){
	df = power.t.test( n=grid$n_subjects[i], delta=grid$delta[i], sig.level = 0.05/n_tests, type="paired")
    data.frame(delta = grid$delta[i], 
    	n_subjects = grid$n_subjects[i], 
    	power= df$power)
	})
res = do.call("rbind", res)

ggplot(res, aes(delta, power, color=factor(n_subjects))) + 
	geom_line(size=1.2) + 
	theme_classic() + 
	theme(aspect.ratio=1) + 
	ylim(0,1) + 
	xlab("Cohen's d") +
	ylab("Power")


# Using observed beta and se_beta
#################################

library(variancePartition)
data(varPartData)

fit = dream( geneExpr, ~ Batch , info)
fit = eBayes(fit)

# Get beta and se_beta for each gene
tab = topTable( fit, coef="Batch2", number=Inf )



# Given the observed effect sizes, t-statistic is t = beta/se.
# Increasing the sample size from n1 to to n2 will decrease 
# the se by a factor of sqrt(n1/n2).
# Compute power assuming effect sizes and 
# estimated standard  errors are accurate
powerWithMoreSamples = function(tab, n1, n2, n_covariates){

	tab$se = with(tab, logFC/t)

	res = lapply(seq(length(n2)), function(i){

		# Degrees of freedom of test
		df = n2[i] - n_covariates

		# compute new t-statistic
		tab$t.new = with(tab, logFC / (se*sqrt(n1/n2[i])))

		# compute new p-value
		tab$P.Value.new = with(tab, 2*pt(abs(t.new), df, lower.tail=FALSE))

		# compute DE genes
		nDE = sum(p.adjust(tab$P.Value.new, "bonferroni") < 0.05)

		data.frame(SampleSize = n2[i], nDE)
	})
	res = do.call(rbind, res)

	res
}


# original sample size
n1 = nrow(info)

# Increased sample size
n2 = seq(100, 1000, by=100)

res = powerWithMoreSamples( tab, n1, n2, 4)

ggplot(res, aes(SampleSize, nDE)) +
	geom_point() +
	geom_line() +
	theme_classic() + 
	theme(aspect.ratio=1) + 
	xlab("Sample size") +
	ylab("Number of differentially expressed genes (expected)") 


# Chinwe's dataset
###################

tab = read.table("~/Downloads/topTableAbs1.poweranalysis5.4.csv", sep=",", header=TRUE, row.names=1)

# original sample size
n1 = 44

# Increased sample size
n2 = c(n1, seq(50, 300, by=20))

res = powerWithMoreSamples( tab, n1, n2, 6)

ggplot(res, aes(SampleSize, nDE)) +
	geom_point() +
	geom_line() +
	theme_classic() + 
	theme(aspect.ratio=1) + 
	xlim(0, max(n2)) + 
	xlab("Sample size") +
	ylab("Number of differentially expressed genes (expected)") 



# Converting between t-statistic and Cohen's d
# Where a and b are the sample sizes in each group
# d = t * sqrt((a+b)/(a*b))
a = 16
b = 22
tab$d = with(tab, t*sqrt((a+b)/(a*b)))






