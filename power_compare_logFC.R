# Gabriel Hoffman
# December 5, 2019

# Compute power for comparing two logFC values
# Motivated by Nadine's synergy analysis
# the sample size enters only through the standard error of the logFC's
# We use the factor that increasing the sample size by a factor of F
#	will decrease the standard error of the logFC by a factor of sqrt(F)
# Were we show results for testing the absolute difference in logFC values
# We also account for multiple testing


#' @param sig1 standard error of estimated for logFC1
#' @param sig2 standard error of estimated for logFC2
#' @param N sample size for which sig1 and sig2 were estimated
#' @param N_other other sample size to evaluate power for
#' @param alpha nominal false positive rate for a single test
#' @param n_tests number of tests
power_compare_logFC = function( sig1, sig2, N, N_other, alpha=0.05, n_tests=20000){

	# absolute difference in logFC
	d = seq(0, 2, length.out=1000)

	# cutoff after multiple testing
	alpha_multiple = alpha / n_tests

	# N_other/N is the relative sample size
	df = lapply( N_other/N, function(n_scale){
		
		# variance of logFC1 - logFC2 is sig1^2 + sig2^2
		# since standard error of the mean is inversely proportional to sqrt(N)
		# 	multiplying the sample size by F decrease the SE by sqrt(F)
		# On the variance scale, this corresponds to dividing by n_scale
		sigSq = (sig1^2 + sig2^2) / n_scale

		# find the cutoff corresponding to the alpha value
		cutoff = qnorm( alpha_multiple/2, 0, sd=sqrt(sigSq), lower.tail=FALSE)

		# power for one sided test, negative
		p1 = pnorm(-1*cutoff, d, sqrt(sigSq))

		# power for second side of test, positive
		p2 = 1-pnorm(cutoff, d, sqrt(sigSq))

		# total power is sum of p1, p2
		data.frame(n_scale, d, power=p1+p2)
	})
	df = do.call("rbind", df)

	# make plot
	ggplot(df, aes(d, power, color=as.factor(n_scale*N))) + geom_line() + theme_bw(14) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + ylim(0, 1) + scale_color_discrete(
		"Samples") + xlab(bquote(abs(logFC[observed] - logFC[expected]))) + ggtitle("Power versus difference in logFC")
}


library(ggplot2)

power_compare_logFC(sig1	= 0.1, 
					sig2	= 0.2, 
					N		= 3, 
					N_other	= c(4, 8, 12, 18, 24), 
					alpha	= 0.05, 
					n_tests	= 20000)







