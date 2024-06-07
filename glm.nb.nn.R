
# Non-negative Negative Binomial regression for deconvolution
# Gabriel Hoffman
# June 7, 2024

library(glmnet)
library(MASS)
library(gamlss)

fit1 <- glm.nb(Days ~ Sex + Age, data = quine)

glm.nb.nn = function(X, y, weights = NULL, offset=NULL){

	# objective function
	ll = function(theta.a){
		# transform theta from a real to positive value
		theta = mgcv::notExp(theta.a)

		fit <- glmnet(X, y, 
			offset = offset,
			weights = weights,
			lower.limits = 0, 
			lambda = 0, 
			family = negative.binomial(theta = theta))
		fit$dev.ratio
	}

	# optimize
	obj = optimize(ll, interval=c(-1e4, 1e4), maximum=TRUE)

	# back-transform positive to real
	theta.hat = mgcv::notExp(obj$maximum)

	# final fit
	fit <- glmnet(X, y, 
		offset = offset,
		weights = weights,
		lower.limits = 0, 
		lambda = 0, 
		family = negative.binomial(theta = theta.hat))

	beta.hat = as.matrix(coef(fit))
	beta.hat = as.vector(beta.hat)
	names(beta.hat) = rownames(coef(fit))

	# remove intercept term
	beta.hat = beta.hat[-1]

	list(beta = beta.hat / sum(beta.hat), theta = theta.hat)
}


X = model.matrix(~ Sex + Age -1 , data = quine)
y = quine$Days

fit = glm.nb.nn(X, y)



library(MuSiC)

GSE50244.bulk.eset = readRDS('~/Downloads/GSE50244bulkeset.rds')
EMTAB.sce = readRDS('~/Downloads/EMTABsce_healthy.rds')
XinT2D.sce = readRDS('~/Downloads/XinT2Dsce.rds')




bulk.mtx = exprs(GSE50244.bulk.eset)

# Estimate cell type proportions
Est.prop.GSE50244 = music_prop(bulk.mtx = bulk.mtx, sc.sce = EMTAB.sce, clusters = 'cellType', samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma','acinar', 'ductal'), verbose = F)
names(Est.prop.GSE50244)

Est.prop.GSE50244$Est.prop.weighted[1:4,]
Est.prop.GSE50244$Est.prop.allgene[1:4,]

library(dreamlet)

EMTAB.sce$static = "a"

pb = aggregateToPseudoBulk(EMTAB.sce, 
					sample_id = "static",
					cluster_id = "cellType")
X_ref = lapply( assayNames(pb), function(x){
		as.matrix(assay(pb, x))
		})
X_ref = do.call("cbind", X_ref)
colnames(X_ref) = assayNames(pb)

sce = SingleCellExperiment(assays = list(counts = X_ref))

X_ref_norm = computeLogCPM(sce)

gene_id = intersect(rownames(X_ref_norm), rownames(bulk.mtx))
X = X_ref_norm[gene_id,c('alpha', 'beta', 'delta', 'gamma','acinar', 'ductal')]
y = bulk.mtx[gene_id,3]

fit = glm.nb.nn(X, y, weights=log(y+2))





