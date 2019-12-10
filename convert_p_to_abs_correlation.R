# Gabriel Hoffman
# December 10, 2019
#
# Given p-value and sample size, get absolute correlation value


# Given a p-value (p) and the sample size (n) compute
# the absolute correlation that gives this p-value
convert_p_to_abs_correlation = function(p, n){
	# Compute t-value
	tstat = qt(p/2, n-2, lower.tail=FALSE)

	# function to optimize
	f = function(r, tstat){
		(tstat - r*sqrt(n-2) / sqrt(1-r^2))^2
	}

	# find correlation (r) that minimizes loss
	fit = optimize(f, c(1e-8, 1-1e-8), tstat=tstat)

	# get correlation value
	fit$minimum
}


# simulation
n = 100
x = rnorm(n)
y = rnorm(n)

cor.test(x,y)

r = cor(x,y)

tstat = r*sqrt(n-2) / sqrt(1-r^2)

pValue = 2*pt(abs(tstat), n-2, lower.tail=FALSE)




# Convert p-value to correlation
convert_p_to_abs_correlation( pValue, n )

# compuare to absolute correlation
abs(r)