# Gabriel Hoffman
# March 31, 2020


library(variancePartition)
library(PRROC)
library(BiocParallel)
# register(SnowParam(4))

# set.seed(1)

n = 100
p = 1000
n_de = 0
n_batch =12

# Generate study design
info = data.frame( Individual = rep(paste0("ID", 1:n), n_batch),
	code = rep(1:n, n_batch), 
	Batch=NA, Disease=NA)

info$Disease[] = "case"
info$Disease[info$code <=n/2] = "control"
for( ID in info$Individual ){
	info$Batch[info$Individual == ID] = paste0("Batch", 1:n_batch)
}
info$Batch = factor(info$Batch)

# randomly drop K samples 
K = n*n_batch*.8
exclude = sample.int(nrow(info), K)
info = info[-exclude,]


# Generate gene expression
Y = lapply( 1:p, function(i){
	if( i < n_de ){
		beta = rnorm(1, 0, 5)
	}else{
		beta = 0
	}
	eta_disease = model.matrix( ~ Disease, info) %*% c(0, beta)
	eta_batch = model.matrix( ~ 0+Batch, info) %*% rnorm(n_batch, 0, 3)
	eta_ID = model.matrix( ~ 0+Individual, info) %*% rnorm(n)
	y = eta_disease + eta_batch + eta_ID + rnorm(nrow(info), 0, .03)

	t(y)
	})
Y = do.call("rbind", Y)




# comment out for now
# # fit variancePartition model
# vp = fitExtractVarPartModel( Y, ~ (1|Disease) + (1|Batch) + (1|Individual), info)

# plotVarPart( vp )

# # fit model
# fit = dream( Y, ~Disease + Batch + (1|Individual), info)
# res = topTable(fit, coef="Diseasecontrol", number=Inf, sort.by="none")

# # show performance
# # This is a step toward validating the model, but here just on fully data
# pr <- pr.curve( -log10(res$P.Value[1:n_de]), -log10(res$P.Value[-c(1:n_de)]), curve=TRUE, rand.compute=TRUE )

# plot(pr, rand.plot=TRUE)





# library(edmcr)

# # Estimate batch effect directly
# ################################

# batchOffsets = lapply( 1:p, function(i){

# 	idx = combn(n_batch, 2)
# 	C = matrix(NA,n_batch,n_batch)
# 	C_weight = C

# 	for(h in 1:ncol(idx) ){
# 		j = paste0("Batch", idx[1,h])
# 		k = paste0("Batch", idx[2,h])
		
# 		# get samples in batches j and k
# 		idx_j = which(info$Batch %in% j)
# 		idx_k = which(info$Batch %in% k)

# 		# only keep individuals present in both batches
# 		keepID = intersect(info$Individual[idx_j], info$Individual[idx_k])
# 		idx_id = which(info$Individual %in% keepID)

# 		# only keep samples in batch j and k where the 
# 		# individuals are present in both batches
# 		isec = intersect(idx_id, c(idx_j, idx_k) )

# 		if( i == 1){
# 			cat(k, j, length(isec), "\n")
# 		}
# 		# if there is overlap
# 		if( length(isec) > 2){
# 			# fit paired model of the different between samples from the
# 			# same individual
# 			fit_local = lm(Y[i,isec] ~ Batch + Individual, droplevels(info[isec,]))
# 			beta = coef(fit_local)[2]
# 			weight = 1/coef(summary(fit_local))[2,2]
# 		}else{
# 			beta = NA
# 			weight = 0
# 		}

# 		C[idx[1,h], idx[2,h]] = beta
# 		C_weight[idx[1,h], idx[2,h]] = weight
# 	}
# 	diag(C) = 0
# 	diag(C_se) = 0
# 	C[lower.tri(C)] = t(C)[lower.tri(C)]
# 	C_weight[lower.tri(C_weight)] = t(C_weight)[lower.tri(C_weight)]

# 	# C[1,2] = C[2,1] = NA
# 	# # matrix completion for distance matrix
# 	# # estimate missing values 
# 	if( any(is.na(C)) ){
# 		# C_est = edmc(C,"dpf", d=1)$D

# 		dc = npf(C, A = C_weight, d=1 )
# 		C_est = dc$D
# 		bv = dc[3][[1]][,1]

# 	}else{
# 		C_est = C
# 	}
# 	# C_est = C

# 	# # project distances onto a line
# 	# # there is an issue with the sign of MDS being inconsistent
# 	batch_values = t(cmdscale( C_est, k=1))
# 	# batch_values = eigen(C_est)$vectors[,1]

# 	fctr = ifelse(sign(C_est[2,1]) == sign(batch_values[2] - batch_values[1]),
# 				1, -1 )

# 	batch_values * fctr
# 	# b = C_est[,1]

# 	# if( cor(a,b) < 0){
# 	# 	stop(i)
# 	# }

# 	})
# batchOffsets = do.call(rbind, batchOffsets)
# colnames(batchOffsets) = levels(info$Batch)




# # Apply the offsets
# ###################

# Y_corrected = matrix(NA, nrow(Y), ncol(Y))

# for( batch in levels(info$Batch) ){

# 	idx = (info$Batch == batch)

# 	# add offset to each sample in each batch
# 	# use different offset (i.e. row) for each gene
# 	Y_corrected[,idx] = Y[,idx] - batchOffsets[,batch]
# }


# Use linear mixed model to analyze the repeated measures
# If we measure an individual mutliple times, compress
# that into a single value per individual
# library(lme4)

# R code to extrat estimate at Individual level
# Y_corrected = lapply( 1:p, function(i){

# 	# Get value for each individual
# 	fit = lmer(Y[i,] ~ Batch + (1|Individual), info)
#  	t(ranef(fit)$Individual)
# })
# Y_corrected = do.call("rbind", Y_corrected)


# Simpler variancePartition call
fitList = fitVarPartModel( Y, ~ Batch + (1|Individual), info, fxn = function(fit){ t(ranef(fit)$Individual) }, showWarnings=FALSE)

Y_corrected = do.call(rbind, fitList)





# estimate batch offset from original data
# lm( t(Y[1,,drop=FALSE]) ~ Batch, info)

# # estimate batch offset from corrected data
# # Observe *no* batch effect
# # here, with an intercept term
# lm( t(Y_corrected[1,,drop=FALSE]) ~ Batch, info)

# # here no intercept terms so values are equal
# lm( t(Y_corrected[1,,drop=FALSE]) ~ 0+Batch, info)


# original values
# Y[1:3, info$Individual == "ID2"]

# # corrected values
# Y[1:3, info$Individual == "ID2"] - batchOffsets[1:3,]

# # Also, corrected values
# # But notice that the values are not exactly the same
# # The method does not set the values from the pair of samples to be equal
# # For a given gene, it computes the deviation between all pairs of individuals, 
# # then creates an offset so that the mean deviation between the samples from the same individuals is zero
# Y_corrected[1:3, info$Individual == "ID2"] 

# Plot coloring by individual
par(mfrow=c(1,2))
# PCA of original data
dcmp = prcomp( t(Y))
(dcmp$sdev^2 / sum(dcmp$sdev^2))[1:2]
plot(dcmp$x[,1:2], col=info$Individual, main="Original data")

# PCA of corrected data
# now disease state stands out
dcmp = prcomp( t(Y_corrected))
(dcmp$sdev^2 / sum(dcmp$sdev^2))[1:2]
plot(dcmp$x[,1:2], col=info$Individual, main="Corrected data")


# vp_orig = fitExtractVarPartModel( Y, ~ (1|Disease) + (1|Batch) + (1|Individual), info)
# vp_corrected = fitExtractVarPartModel( Y_corrected, ~ (1|Disease) + (1|Batch) + (1|Individual), info)

# plot(vp_orig$Batch, vp_corrected$Batch)
# abline(0,1, col="red")




# # NOTE
# Here I used the ideal example of complete observations in all datsets.  
# In this case all distances between batches are guaranteed to sum (i.e. d(1,3) = d(1,2) + d(2,3)) and all pairs are observed.
# Also in this case the observe batch effect from the corrected values is guaranteed to be zero.
# I started with this to make sure it works in the best case scenerio.
# Your dataset corresponds to the case where some samples are dropped from Y_corrected.  In that case, the observed batch effect of the corrected data will not be exactly zero, because theere is incomplete overlap between the "training" and "testing" sets.




# For gene i, compare batch j and batch k on the samples from paired individuals that overlap
# Create a matrix C where C[j,k] is the batch effect between j and k. Set diagonal values to zero.

# 1) Fit a linear model on the subset of data to estimate C[j,k]
# 	fit = lm(y ~ Batch)
# 	Repeat for all pairs of batches

# 2) If some entries in C are not observed (i.e. NA) because no paired samples overlap between batches j and k, use Euclidean Distance Matrix Completion (edmcr package) to fill in these values.

# 	# matrix completion for distance matrix
# 	# this assumes that batches lay on a one dimension line
# 	# so that the distance from batch 1 to 3 is the sum of 1,2 + 1,3
# 	# When the data is completely observed and all pairs overlap in all batches,
# 	# this assumption is satisfied *exactly*.
# 	# With incomplete overlap, this assumption can be used to reconstruct missing distances
# 	# If all pairs are observed, then skip this
# 	library(edmcr)
# 	C_est = edmc(C,"dpf", d=1)$D

# 	# Given this complete pairwise distance matrix, 
# 	# create a vector batch_values where batch_values[h] 
# 	# should be added to samples in batch h to remove the batch effect.
# 	# As above, this assumes that the batches effects lay on a 1 dimensional line
# 	# When the data is completely observed and all pairs overlap in all batches, this will give the offset terms that you expect intuitively
# 	batch_values = cmdscale( C_est, k=1)
# 	batch_values = batch_values - min(batch_values) # so min offset is zero

# 3) Apply the batch effect correction by adding batch_values[h] to the expression values from batch h











