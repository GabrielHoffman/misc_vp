

1) Power 
	150 train, 40 test

	power for prediction
	assume effect sizes

2) Power for causal networks: Noam
	



3) 

Cohens d is the effect size divided by its standard error

library(ggplot2)
library(pwr)

# correct for 15K genes for multiple testing
n_tests = 15000
n_subjects = 500


res1 = lapply( seq(0, 1, length.out=1000), function(delta){
    res = power.t.test( n=n_subjects, delta=delta, sig.level = 0.05/n_tests)
    data.frame(power= res$power, delta, Test = "First set of assumptions"    )
    })
res1 = do.call("rbind", res1)

ggplot(res1, aes(delta, power)) + geom_line() + theme_bw() + theme(aspect.ratio=1) + ylim(0,1) + xlab("Cohen's d")


Cohens d is the effect size scaled by the within group standard deviation


library(Rfast)
library(ggplot2)
library(glmnet)
library(parallel)

# size of training set
n1 = 5000
# size of testing set
n2 = 100
# number of features
p = 100

# simulate over a range of fraction of variance explaind by predictions
df = mclapply( seq(0, .95, length.out=10), function(FVE){

	df = lapply(c(100, 200, 300, 500, 1000), function(n1){

		# simulate features
		X = matrnorm(n1+n2, p)

		df = lapply(1:10, function(f){

			# simulate coefficients for each feature
			beta = rnorm(p)

			# create response value with no noise
			eta = X %*% beta

			# add noise for a given FVE value
			# so that squared correlation between true and observed value is FVE
			# cor(y, eta)^2
			if( FVE == 0){
				y = rnorm(n1+n2)
			}else{
				y = scale(eta) + rnorm(n1+n2, 0, sd=sqrt(1/FVE - 1))
			}
			y = scale(y)
			# y = rbinom(length(eta), 1, prob=inv.logit(eta))		
		
			# fit = cv.glmnet(X[1:n1,], y[1:n1], family="binomial")

			# for lasso regression with cross validation
			fit = cv.glmnet(X[1:n1,], y[1:n1])

			# get the predicted values for the testing set
			# based on the coefficents from the best model 
			y_hat = predict(fit, newx = X[(n1+1):(n1+n2),], s = "lambda.1se")

			# compute the root mean squared error for preditions
			rMSE = sqrt(mean((y_hat - y[(n1+1):(n1+n2)])^2))

			data.frame(n1, n2, FVE, rMSE)
		})
		df = do.call(rbind, df)

		# report the mean and standard deviation of rMSE estimates over 5 simulations
		data.frame(n1 = df$n1[1], rMSE = mean(df$rMSE), se = sd(df$rMSE), FVE=FVE)
	})
	do.call(rbind, df)

}, mc.cores=10)
df = do.call(rbind, df)

df$n1 = factor(df$n1, sort(unique(df$n1)))


ggplot(df[df$n1=="500",], aes(FVE, rMSE, color=n1)) + theme_bw() + theme(aspect.ratio=1) + xlab("Fraction of variance explained by predictors") + geom_errorbar(aes(x=FVE, ymin=rMSE - se, ymax=rMSE+se), width=0) + geom_point() + xlim(0, 1)  + geom_smooth(method="loess", se=FALSE)
















