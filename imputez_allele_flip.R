
# Matching alleles between GWAS and reference panel
suppressPackageStartupMessages({
library(VariantAnnotation)
library(tidyverse)
library(imputez)
library(progress)
library(Rfast)
})

# Read GWAS results
fileGWAS = "/sc/arion/projects/roussp01b/resources/databases/gwas/sz3/PGC3-SCZ.txt.gz"
df_GWASIn = read_tsv(fileGWAS, show_col_types=FALSE)

# filter out NA p-values
# convert p-values and ORs to z-scores
df_GWASIn = df_GWASIn %>% 
		filter(!is.na(P), INFO > 0.8) %>%
		mutate(z.stat = qnorm(P / 2,lower.tail = FALSE) * sign(OR - 1)) %>%
		dplyr::select(SNP, A1, A2, OR, P, z.stat) 

# Read VCF
file = "/sc/arion/projects/CommonMind/christian/CMC-VNTR/files/imputez_chromwise/chr_22/CMC_imputed_hg38_maf01_DLPFC_reheader_qc4_22.vcf"
vcfIn <- readVcf(file, "hg38")

# keep only SNV/INDEL siets with one ALT allel
keep = isSNV(vcfIn, singleAltOnly=TRUE) | isIndel(vcfIn, singleAltOnly=TRUE)
vcfIn = vcfIn[keep,]

# assumes only biallelic sites
df_position = tibble(SNP = rownames(vcfIn), 
				REF = as.character(ref(vcfIn)),
				ALT = as.character(unlist(alt(vcfIn))))

# merge GWAS results and VCF alleles
df_GWAS = inner_join(df_GWASIn, df_position, by="SNP") 

# only keep SNPs with same allelic set, regardless of order
df_GWAS = df_GWAS %>%
	filter(all(sort(c(A1, A2)) == sort(c(REF, ALT))	))

# compute signFlip based on matching A1 to REF
# then apply flip
df_GWAS = df_GWAS %>%
	mutate(signFlip = ifelse(A1 == REF, 1, -1)) %>%
	mutate(z.stat.final = z.stat * signFlip)

# Get X_geno as dosage
X_geno.tmp = geno(vcfIn)$GT
X_geno = matrix(0, nrow=nrow(vcfIn), ncol=ncol(vcfIn), dimnames = dimnames(geno(vcfIn)$GT))
X_geno[X_geno.tmp=="0/0"] = 0
X_geno[X_geno.tmp=="0/1"] = 1
X_geno[X_geno.tmp=="1/0"] = 1
X_geno[X_geno.tmp=="1/1"] = 2


run_impute_geno = function(z, X, query, lambda = NULL, windowSize=200, quiet=FALSE){

	if (!quiet) {
		pb <- progress_bar$new(
		format = "  imputing [:bar] :percent eta: :eta",
		total = length(query), clear = FALSE, width = 60)
	}

	df = lapply( query, function(id){

		# get SNP index and window
		i = match(id, names(z))
		idx = seq(max(1, i-windowSize), min(nrow(df_GWAS), i + windowSize))

		# compute corrleation
		X_sub = t(X[names(z)[idx],])
		Sigma = cora( X_sub )

		if( is.null(lambda) ){
			lambda.est = corpcor::estimate.lambda(scale(X_sub), verbose=FALSE)
		}else{
			lambda.est = lambda
		}

		# get location of target SNP
		j = match(id, names(z)[idx])

		# Impute
		res = imputez( z[idx], Sigma, j, lambda=lambda.est)
		res$z.observed = z[idx[j]]

		if (!quiet) pb$tick()

		tibble(res)
		})
	if (!quiet) pb$terminate()
	bind_rows(df)
}



z = df_GWAS$z.stat.final
names(z) = df_GWAS$SNP


windowSize = 300

idx = sample.int(nrow(df_GWAS), 1000)
query = df_GWAS$SNP[idx]

df = run_impute_geno(z, X_geno, query, lambda = .1, windowSize = windowSize)

# cor all 
cor(df$z.observed, df$z.stat)

# corr good quality
with(df %>% filter(r2.pred > 0.8), cor(z.observed, z.stat))

limits = max(abs(c(df$z.observed, df$z.stat)))

ggplot(df, aes(z.observed, z.stat, color=r2.pred)) +
	geom_point() +
	theme_classic() +
	theme(aspect.ratio=1) +
	scale_color_gradient(high="black", low="red", limits=c(0,1)) +
	geom_abline(intercept=0, slope=1) +
	scale_x_continuous(limits=c(-limits, limits)) +
	scale_y_continuous(limits=c(-limits, limits)) +
	xlab("Observed z-statistic") +
	ylab("Imputed z-statistic")







## OLD CODE
###########





to_str = function(allele){
	unlist(allele, recursive=TRUE)
	# sapply(unlist(allele), function(x){
	# 		x = as.character(x)
	# 	ifelse( length(x) == 1, x, NA)
	# 	})
}

# extract ALT allele
df_alt_allele = rowRanges(vcfIn)$ALT %>% 
				unlist %>%
				data.frame(ID = rownames(vcfIn), ALT=.) %>%




df_position = data.frame(ID = rownames(vcfIn), 
				REF = as.character(ref(vcfIn)),
				ALT = as.character(unlist(alt(vcfIn))))




rowRanges(vcfIn) %>%
	data.frame(ID = rownames(.))

	head(1000) %>%
	as.data.frame %>%
	mutate(ALT1 = unlist(ALT, recursive=TRUE)) %>%
	tibble




df = data.frame(REF = unlist(rowRanges(vcfIn)$REF), ALT = unlist(rowRanges(vcfIn)$ALT))

alt_allele = lapply(rowRanges(vcfIn)$ALT, function(x){
	x = as.character(x)
	data.frame(ALT = x[1], n_alleles = length(x))})
alt_allele = do.call(rbind, alt_allele)



df_position = rowRanges(vcfIn[1:10,]) %>% 
					as.data.frame %>%
					tibble(ID = rownames(.), .) %>%
					dplyr::select(ID, CHROM = seqnames, REF, as.character(ALT))

df_GWAS = df_GWASIn %>%
		filter(SNP %in% rownames(vcfIn))


df_merge














library(data.table, quietly = TRUE)

library(dplyr, quietly = TRUE)

library(tidyr, quietly = TRUE)

library(ggplot2, quietly =TRUE)

library(stringr, quietly = TRUE)

library(ggpubr, quietly = TRUE)

#

file = "/sc/arion/projects/roussp01b/resources/databases/gwas/sz3/PGC3-SCZ.txt.gz"

dosage = "/sc/arion/projects/CommonMind/christian/CMC-VNTR/files/imputez_chromwise/chr_22/CMC_imputed_hg38_maf01_DLPFC_reheader_qc4_22.vcf"

allele_check = as.data.frame(fread(dosage))

df_z_obs = as.data.frame(fread(file,fill = TRUE))

df_z_obs$P[df_z_obs$P == 1.0000] <- 0.9999 

#

df_z_obs=merge(df_z_obs,allele_check,by.x='SNP',by.y='ID') # Allele Check

df_z_obs$OR2<-ifelse(df_z_obs$A1==df_z_obs$REF,df_z_obs$OR,1/df_z_obs$OR) # Match the alleles

df_z_obs$Z<-qnorm(df_z_obs$P,lower.tail = FALSE)

z=df_z_obs[df_z_obs$CHR==chr,]

z$Z_2<-ifelse(sign(z$OR2 >= 1.0)=="0",z$Z*-1,z$Z*1) #if OR not greater than 1, multiply Z-score by -1

###

df_testing_500kb<-readRDS("/sc/arion/projects/CommonMind/christian/imputez_snps_for_Gabe.RDS")

df2<-merge(df_testing_500kb,z,by.x='ID',by.y='SNP')

df2$pval.stat<-2*pnorm(-abs(df2$z.stat),lower.tail = TRUE)

df2<-df2[df2$REGION=='DLPFC',]

r = with(df2, cor(Z_2, z.stat))

### Plot ###

ggplot(df2, aes(Z_2, z.stat, color=r2.pred)) +
  
  geom_point() +
  
  theme_classic() +
  
  theme(aspect.ratio=1) +
  
  geom_abline(color="red") + 
  
  xlab("Observed z") +
  
  ylab("Imputed z") + 
  
  labs(title="Compare observed and imputed SNPs",subtitle = "Window Size 500kb") +
  
  scale_color_gradient(low="lightblue", high="black", limits=c(0,1)) +
  
  annotate(geom="text", x=-2, y=2, label=paste0("R=",format(r, digits=4)))