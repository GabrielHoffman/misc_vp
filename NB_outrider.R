
# Gabriel Hoffman
# Sept 15, 2021

library(MASS)
library(limma)
library(variancePartition)
library(edgeR)
library(parallel)

n_genes = 20000
n_samples = 100

mu = 10^runif(n_genes, -1, 4)

hist(mu)

geneCounts = lapply(1:n_samples, function(i){
	sapply(mu, function(x) rnegbin(1, x, theta = 1e3))
})

geneCounts = do.call(cbind, geneCounts)
colnames(geneCounts) = paste0('S', 1:n_samples)
rownames(geneCounts) = paste0('G', 1:n_genes)

colSums(geneCounts)


dge = DGEList(counts = geneCounts)
dge = calcNormFactors(dge)
keep = filterByExpr(dge)
dge = dge[keep,]

vobj = voom( dge, model.matrix(~1, info),plot=TRUE)


info = data.frame(Donor = colnames(geneCounts))
formula = ~1

# extrema = function(vobj, formula, info){

# Standard
fit = dream(vobj, formula, info)

resid = residuals(fit)

zScores = apply(resid, 1, function(x){
	scale(x)
	})

hist(zScores)

pv = 2*pnorm(abs(zScores), lower.tail=FALSE)
z.fdr = p.adjust(pv, "fdr")
sum(z.fdr < 0.05)


fit_robust = lmFit(vobj, model.matrix(formula, info), method="robust")

resid = residuals(fit_robust, vobj)

zScores = apply(resid, 1, function(x){
	scale(x)
	})

hist(zScores)

pv = 2*pnorm(abs(zScores), lower.tail=FALSE)
z.fdr = p.adjust(pv, "fdr")
sum(z.fdr < 0.05)


tabList = mclapply(unique(info$Donor), function(id){

	message("\r", id, "    ", appendLF=FALSE)

	dsgn = model.matrix(~ I(Donor==id), info)

	fit = lmFit(vobj, dsgn, method="robust")
	fit = eBayes(fit, robust=TRUE)

	topTable(fit, coef=2, number=Inf, sort.by="none")

}, mc.cores=4)


pv = sapply(tabList, function(tab) tab$P.Value)
z.fdr = p.adjust(pv, "fdr")
sum(z.fdr < 0.05)

tab = do.call(rbind, tabList)

 with(tab, plot(logFC, -log10(P.Value)))

  with(tab, plot(AveExpr, logFC))







library(OUTRIDER)

ods <- OutriderDataSet(countData=geneCounts)

# filter out non expressed genes
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)

# run full OUTRIDER pipeline (control, fit model, calculate P-values)
ods <- OUTRIDER(ods)








