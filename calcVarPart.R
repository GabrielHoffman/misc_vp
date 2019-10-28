# Gabriel Hoffman
# October 28, 2019
# Variance partition fraction for GLM with linear, logistic or probit regression
# Reports variance fraction on the linear scale (i.e. eta)
# rather than the observed scale (i.e. y)
# For linear regression with an identity link, these are the same
# but for logit/probit, the fractions are not well defined on the observed scale
#	[See development code for attempts]
# The residual variance is the variance of the distribution function of the link
# 	logit -> logistic distribution -> variance is pi^2/3
#	probit -> standard normal distribution -> variance is 1
#
# Proposed by
# McKelvey and Zavoina. A statistical model for the analysis of ordinal level dependent variables. The Journal of Mathematical Sociology 4(1) 103-120 https://doi.org/10.1080/0022250X.1975.9989847

# Also see
# DeMaris. Explained Variance in Logistic Regression: A Monte Carlo Study of Proposed Measures. Sociological Methods & Research 2002 https://doi.org/10.1177/0049124102031001002


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,glm-method
setMethod("calcVarPart", "glm",
function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...){

	N = nrow(fit$model)

	# Compute eta for each term
	# predicted value in linear space for each term
	Eta = predict(fit, type="terms")

	# compute residual term for each link
	distVar = switch( fit$family$link , 
			"logit" = (pi^2)/3,
			"probit" = 1,
			"identity" = var(fit$residuals) )

	if( is.null(distVar) ){
		stop("glm link not supported: ", fit$family$link )
	}

	# variance on linear scale
	# get variance of each term
	# append with variance due to link function
	var_term = c(apply(Eta, 2, var), distVar)

	# compute fraction
	varFrac_linear = var_term / sum( var_term )
	names(varFrac_linear) = c(colnames(Eta), "Residuals")

	return( varFrac_linear )
})



# CODE FROM development

	# g^{-1}( mu + eta_j) / var(Y)

	# compute denominator
	#####################

	# # get response
	# y = as.numeric(fit$model[,attr(fit$terms,"response")])

	# # denominator
	# var_Y = var( y )

	# variance on response scale
	# g_inv = fit$family$linkinv

	# Y_hat = apply( Eta, 2, function(x){
	# 	g_inv( x + mu )
	# })
	# # C = cov(Y_hat)
	# # sum(eigen(C)$values)

	# v = apply(Y_hat, 2, var)


	# v / var_Y

	# denom = ifelse( var_Y > sum(v), var_Y, sum(v))



	# varFrac_response = apply(Y_hat, 2, var) / var_Y


	# apply(Y_hat, 2, var) / var_Y


	# # deviance 
	# deviance = anova(fit)[['Resid. Dev']]
	# 1 - (fit$deviance / deviance[1])





	# baseline
	# mu = attr(Eta,"constant")

	# par(mfrow=c(1,3))
	# plot(y, rowSums(Eta) + mu)
	# plot(y, predict(fit))
	# plot(y, predict(fit, type="response"))

	# attr(fit$terms,"term.labels")



# dcmp = svd(scale(Y_hat))

# A = sweep(dcmp$u, 2, dcmp$d^-1,FUN="*")

# cor(t(A) %*% Y_hat)




# dcmp = svd(scale(Y_hat))
# A = sweep(dcmp$v, 2, dcmp$d^-1,FUN="*")
# cor(scale(Y_hat) %*% (A))




# dcmp = svd(scale(Y_hat, scale=FALSE))
# A = sweep(dcmp$v, 2, dcmp$d^-1,FUN="*")
# cor(scale(Y_hat, scale=FALSE) %*% A)
# cv = colVars(scale(Y_hat, scale=FALSE)%*% A)
# cv
