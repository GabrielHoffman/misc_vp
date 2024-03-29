#! /usr/bin/env Rscript
#
# Gabriel Hoffman
# May 20, 2021
#
# Genotype concordance like NGScheckmate, but using genotype dosage
# instead of read depth.  This code works with VCF generated by array and 
# sequencing, while NGScheckmate fails with genotype data. 
# Includes single or multi-sample VCFs.


# check_vcf_concordance.R --snpBed /sc/arion/projects/H_PBG/REFERENCES/GRCh38/NGSCheckmate/hglft_genome_6702d_cdce40.bed \
# --vcfList /sc/arion/projects/CommonMind/shan/MOLECULAR_PROFILING/ngscheckmate/final_vcfs.lst \
# --nthreads 7 \
# --outfile correlation.tsv.gz


library(getopt)

spec = matrix(c(
      'snpBed', 	'd', 1, "character",
      'vcfList',	'v', 1, "character",
      'nthreads',	'n', 1, "numeric",
      'outfile',	'o', 1, "character"
    ), byrow=TRUE, ncol=4)
opt = getopt(spec)

if( ! file.exists(opt$snpBed) ){
	stop("File does not exist:", opt$snpBed)
}

if( ! file.exists(opt$vcfList) ){
	stop("File does not exist:", opt$vcfList)
}

if( is.null(opt$nthreads)){
	opt$nthreads = 1
}

if( is.null(opt$outfile)){
	stop("Must specify --outfile")
}

message("Loading libraries...")
suppressPackageStartupMessages({
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(parallel)
library(comperes)
library(R.utils)
})

files = read.table(opt$vcfList)$V1

# read in target SNP data
gr = with(fread(opt$snpBed), GRanges(V1, IRanges(V2, V3, name=V4, REF=V5, ALT=V6)))
gr = keepSeqlevels(gr, paste0("chr", 1:22), pruning.mode="coarse")
gr$ID = with(gr, paste0(seqnames, ':', end, ':', REF, ':', ALT))

# Compute dosage from genotype probabilities
getDosage = function(genProb){
	value = sapply(genProb, function(x){
		idx = is.na(x)
		if( length(idx) > 0){
			x[idx] = 0
		}
		# convert probability to dosage
		x[2] + 2*x[3]
		}  )
	M = matrix(value, nrow=nrow(genProb), byrow=FALSE)
	rownames(M) = rownames(genProb)
	colnames(M) = colnames(genProb)
	M
}

# check if file is a bgzipped and tabix
isTabixed = function(file){

	tryCatch({
		h = headerTabix( file )
		TRUE
	}, error = function(e){
		FALSE
	})
}

isSmall = function(file, cutoff = 200){
	# get size of file in Mb
	sizeMB = file.info(file)$size / 1000000
	sizeMB < cutoff
}

# check that files exist
res = sapply( files, function(file){
	if( ! file.exists(file) ){
		stop("VCF file does not exist:\n", file)
	}
	})


# for each VCF
message("Reading VCFs...")
dosageData = mclapply( files, function(file){

	# read VCF
	if( isTabixed(file) & ! isSmall(file) ){
		res = readVcf( file, genome = "GRCh38", param = gr )
	}else{		
		res = readVcf( file, genome = "GRCh38" )
	}

	# read variant locations
	gr2 = rowRanges(res)

	# get Dosage for each variant
	if( 'DS' %in% names(geno(res)) ){
		Dosage = geno(res)$DS

	}else if( "PL" %in% names(geno(res)) ){
		genProb = PLtoGP(geno(res)$PL)
		# genProb = matrix(genProb, nrow=nrow(geno(res)$PL))
		Dosage = getDosage(genProb)

		names(gr2) = gsub('_', ':', names(gr2))
		names(gr2) = gsub('/', ':', names(gr2))

	}else if( "GL" %in% names(geno(res)) ){
		genProb = GLtoGP(geno(res)$GL)
		Dosage = getDosage(genProb)

		names(gr2) = gsub('_', ':', names(gr2))
		names(gr2) = gsub('/', ':', names(gr2))
	}
	rm(res)

	# keep only requested variants
	idx = which(names(gr2) %in% gr$ID)

	# return variant ID's and dosages
	data.frame( ID = names(gr2)[idx], Dosage[idx,,drop=FALSE])
	}, mc.cores=opt$nthreads)

# Merge, keeping all variates
message("Merging dosages...")
df_merge = Reduce(function(x,y) merge(x = x, y = y, by = "ID", all=TRUE), dosageData)


# Compute pairwise correlation
# uses fast method to evaluate all pairs
message("Evaluating concordance...")
C = cor(df_merge[,-1], use="pairwise.complete")

# Convert to pairwise tall matrix
df_pw = mat_to_long(C, "Sample_1", "Sample_2", "correlation")
df_pw = df_pw[!(df_pw[,1] == df_pw[,2]),]
df_pw$Match = df_pw$correlation > 0.9

# count number of non-NA's between each pair of samples
C_counts = crossprod(!is.na(df_merge[,-1]))
df_pw_counts = mat_to_long(C_counts, "Sample_1", "Sample_2", "N")

# Merge correlation and counts
df_pw = merge(df_pw, df_pw_counts, by=c("Sample_1", "Sample_2"))

# write to file
message("Writing to file...")
write.table(df_pw, file=gzfile(opt$outfile), row.names=FALSE, quote=FALSE, sep="\t")

















