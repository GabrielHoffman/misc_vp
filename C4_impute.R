# Gabriel Hoffman
#
# eEpt 14, 2023
# 
# Impute MHC C4 z-statistics using Gaussian approximation


# MHC C4
# https://mccarrolllab.org/resources/resources-for-c4/
 wget --no-check-certificate https://personal.broadinstitute.org/giulio/panels/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz
 ml tabix
 tabix -p vcf MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz


library(VariantAnnotation)
library(Rfast)
library(ggplot2)

file = "/sc/arion/projects/CommonMind/hoffman/ldref/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh38.vcf.gz"
vcf <- readVcf(file, "hg38")

# GT = readGT(file)

# res = t(apply(GT, 1, function(x) as.numeric(as.factor(x))))
# colnames(res) = colnames(GT)


# Sigma = cora(t(res))

# # dcmp = eigen(Sigma)


# lambda = 0.01
# Sigma.shrink = (1-lambda) * Sigma + diag(lambda, nrow(Sigma))

# i = which(rownames(Sigma.shrink) == "C4")
# W = solve(Sigma.shrink[-i,-i], Sigma.shrink[-i,i,drop=FALSE])

# sigSq = Sigma.shrink[i,i] - (crossprod(W,Sigma.shrink[-i, -i]) %*% W)

# 1 - as.numeric(sigSq)



# library(decorrelate)

# ecl = eclairs(t(res[seq(5593-500, 5593+500),]), compute="correlation")



# rowRanges(vcf)["C4",]

# geno(vcf)$GT["C4",]

# get allele codes
alleles = rowRanges(vcf)["C4",]$ALT[[1]]

# convert allele numbers to text codes as factors
GT.C4 = t(sapply(geno(vcf)$GT["C4",], function(x){
	i = as.numeric(unlist(strsplit(x, "\\|")))
	alleles[i]
	}))
df_c4 = data.frame(GT.C4)
df_c4$X1 = factor(df_c4$X1, alleles)
df_c4$X2 = factor(df_c4$X2, alleles)

# sort rows by allele factor
idx = sapply(seq(nrow(df_c4)), function(i) is.unsorted(as.numeric(df_c4[i,])) )

tmp = df_c4$X1[idx]
df_c4$X1[idx] = df_c4$X2[idx]
df_c4$X2[idx] = tmp


X1 = model.matrix(~ 0 + X1, droplevels(df_c4))
X2 = model.matrix(~ 0 + X2, droplevels(df_c4))
X_alleles = cbind(X1, X2)

# NEED to fix 0|1 <=> 1|0 HET coding issue
# imputez does not consider phazing

i = which(rownames(vcf) == "C4")
X_ref = t(apply(geno(vcf)$GT[-i,], 1, function(x){
	y = rep(NA, length(x))
	y[x=="0|0"] = 0
	y[x=="0|1"] = 1
	y[x=="1|0"] = 1
	y[x=="1|1"] = 2
	y
	# as.numeric(as.factor(y))
	}))
colnames(X_ref) = colnames(vcf)

# size of window
i = which(rownames(X_ref) == "C4")
X_ref = X_ref[seq(i-3000,i+2000),]

# append C4 alleles
X_ref_expand = rbind(X_ref, t(X_alleles))



get_R2 = function(Sigma, target, lambda = 0.1 ){

	Sigma.shrink = (1-lambda) * Sigma + diag(lambda, nrow(Sigma))
	# kappa(Sigma)
	# kappa(Sigma.shrink)

	i = which(rownames(Sigma.shrink) == target)
	W = solve(Sigma.shrink[-i,-i], Sigma.shrink[-i,i,drop=FALSE])

	sigSq = Sigma.shrink[i,i] - (crossprod(W,Sigma.shrink[-i, -i]) %*% W)

	1 - as.numeric(sigSq)
}


# get_R2(Sigma, "C4")

library(RhpcBLASctl)
omp_set_num_threads(4)

Sigma.all = cora(t(X_ref_expand))

df_quality = lapply(seq(ncol(X_alleles)), function(i){

	message(i)
	id = colnames(X_alleles)[i]
	keep = setdiff(rownames(X_ref), "C4")
	keep = c(keep, id)

	S = Sigma.all[keep,keep]

	# imputation quality with Gaussian approx
	r2 = get_R2(S, id)

	# max r2 from single marker
	r2.max = max(S[id,rownames(S)!=id]^2)

	data.frame(id = id, 
		r2.max = r2.max, 
		r2 = r2, 
		freq = sum(X_alleles[,i]) / nrow(X_alleles))
})
df_quality = do.call(rbind, df_quality)

ggplot(df_quality, aes(r2.max, r2)) +
	geom_point() +
	theme_classic() +
	theme(aspect.ratio=1) +
	xlim(0, 1) +
	ylim(0, 1) +
	geom_abline(intercept=0, slope=1, linetype="dashed")





S = cora(t(X_ref))
heatmap(S[colnames(X_alleles),]^2, zlim=c(0,1))


S %>%






# SCZ GWAS
##########

# ImputeZ

library(tidyverse)
library(stringr)
library(imputez)

file = "/sc/arion/projects/va-biobank/resources/GWAS/GWAS_INBOX/PGC3plusMVP_CSP572/meta.csp572_mvp006_eur-SCZ_clean.pgc3_scz.rsids.txt.gz"

df_SCZ = read_delim(file, delim=" ", show_col_types = FALSE) %>%
	separate_wider_delim(SNP, names=c("CHR", 'POS', 'A1', 'A2'), ":") %>%
	mutate(POS = as.numeric(POS)) %>%
	arrange(CHR, POS)


df_SCZ_sub = df_SCZ %>%
				filter(MarkerName %in% names(vcf)) %>%
				mutate(z = Effect / StdErr)




# append C4 alleles
X_ref_expand = rbind(X_ref, t(X_alleles))

Sigma.all = cora(t(X_ref))


table(df_SCZ_sub$MarkerName %in% rownames(Sigma.all))

# Append markers without Z-scores
i = rownames(Sigma.all) %in% df_SCZ_sub$MarkerName


df_SCZ_sub2 = as_tibble(rowRanges(vcf) ) %>%
	mutate(SNPID = rownames(vcf)) %>%
	filter(SNPID %in% rownames(Sigma.all)[!i]) %>%
	dplyr::select(seqnames, start, REF, ALT) %>%
	mutate(CHR = gsub("^chr", "", seqnames),
		POS = as.numeric(start), A1 = REF) %>%
	bind_rows(df_SCZ_sub, .)


# Need align alleles


# df_z = run_imputez(df_SCZ_sub$z[1:100], Sigma.all, df_SCZ_sub$MarkerName[1:100])

df_res = lapply(seq(nrow(df_SCZ_sub2)), function(i){
	message(i)

	ids = c(df_SCZ_sub2$MarkerName[i], df_SCZ_sub2$MarkerName[!is.na(df_SCZ_sub2$z)])
	ids = ids[!duplicated(ids)]

	df_z_local = df_SCZ_sub2 %>% 
					filter(MarkerName %in% ids) %>%
					dplyr::select(MarkerName, z)

	z_local = df_z_local$z
	names(z_local) = df_z_local$MarkerName

	Sig.local = Sigma.all[ids,ids]

	res = imputez(z_local, Sig.local, 1, lambda = 0.1)

	res
})
df_res = do.call(rbind, df_res)


df_SCZ_sub2[i,] %>%
	dplyr::select(MarkerName, z)


# Faster using eclairs
######################

library(decorrelate)

X = t(X_ref)[,1:50]

z = rnorm(ncol(X))
names(z) = colnames(X)

imputez(z, cora(X), 1, lambda=0.1)

imputez.eclairs(z, X, 1, lambda=0.1)






	crossprod(W, decorrelate(W, ecl, transpose=TRUE, alpha=1))

	crossprod(W,getCor(ecl) %*% W)



	# crossprod(W,cora(X[,-i]) %*% W)
	(crossprod(W,Sigma.shrink[-i, -i]) %*% W)

	ecl = eclairs(X, compute="correlation", lambda = ecl$lambda)
	Sigma.shrink = getCor(ecl, lambda = lambda)
	cor(X[,1:3])
	Sigma.shrink[1:3, 1:3]
	Sigma.shrink = cora(X)

	W2 = solve(Sigma.shrink[-i,-i], Sigma.shrink[-i,i,drop=FALSE])

	# imputed z-scores
	z_i = crossprod(W2, z[-i])

	# compute standard error for each imputed z-score
	sigSq = Sigma.shrink[i,i] - (crossprod(W2,Sigma.shrink[-i, -i]) %*% W2)
	
	data.frame(	ID = names(z)[i], 
				z.stat = as.numeric(z_i), 
				sigSq = as.numeric(sigSq),
				r2.pred = 1 - as.numeric(sigSq))
}




