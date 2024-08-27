

# cd ~/build2
# cd variancePartition
# # 1.31.16
# git checkout 7074c039dd0b85acb0c6dde53ce26f9f7876089f
# R CMD INSTALL -l /hpc/users/hoffmg01/.Rlib/test .


# cd ~/build2
# cd dreamlet
# # 1.31.16
# git checkout 1c4b727e06a81e0c67e7c76f814b15861e89f6e9
# R CMD INSTALL -l /hpc/users/hoffmg01/.Rlib/test .



library(variancePartition, lib.loc='/hpc/users/hoffmg01/.Rlib/test')
library(dreamlet, lib.loc='/hpc/users/hoffmg01/.Rlib/test')
library(muscat)
library(ExperimentHub)
library(zenith)
library(scater)
library(zellkonverter)
library(dplyr)
library(edgeR)
library(reticulate)

setwd("/sc/arion/scratch/hoffmg01/milos/old")

file = "/sc/arion/projects/CommonMind/mpjanic/CHPLEX/final_csv/final_dreamlet_crumblr/final17.h5ad"

matrx= readH5AD(file, raw=TRUE, use_hdf5=TRUE, verbose=TRUE)
raw <- assay(altExp(matrx, "raw", withDimnames=TRUE, withColData=TRUE), "X")
sce <- SingleCellExperiment(list(counts=raw),
                                colData=colData(matrx),
                                rowData=rowData(matrx),
                                reducedDims=reducedDims(matrx))


pb <- aggregateToPseudoBulk(sce, assay ="counts", cluster_id = "subtype", sample_id="SubID_vS", verbose=FALSE)

saveRDS(pb, file="pb.RDS")


form =  ~ cerad_4_12_24 + scale(age_4_12_24) #+ (1|gender_4_12_24) + (1|set) + log(n_genes)
res.proc = processAssays( pb, form)#, assays="ASTR" )

saveRDS(res.proc, file = "cerad_resproc_subtype_dx.RDS")

contrasts = c(Diff = 'cerad_4_12_24AD - cerad_4_12_24CTR')

form = ~ scale(age_4_12_24) + cerad_4_12_24 + 0  + log(n_genes) #+ (1|set) + (1|gender_4_12_24)
res.dl = dreamlet( res.proc, form, contrasts=contrasts)


saveRDS(res.dl, file = "cerad_resdl_subtype_dx.RDS")


tab = topTable(res.dl, coef="Diff", number=Inf)

# compute FDR within each cell type
res = tab %>%
    as_tibble %>%
    group_by(assay) %>%
    mutate(FDR.within = p.adjust(P.Value, "fdr"))

# counts results
table(res$adj.P.Val < 0.05)
table(res$FDR.within < 0.05)


write.table(res, file="cerad_subtype_dx.csv", sep=",", quote=FALSE)













go.gs = get_GeneOntology(to="SYMBOL")
res_zenith = zenith_gsa(res.dl, coef = 'Diff', go.gs)
write.table(res_zenith, file="cerad_zenith_subtype_dx.csv", sep=",", row.names = FALSE, quote=FALSE)




res.proc = processAssays( pb,  ~ braak_4_12_24 + scale(age_4_12_24) + (1|gender_4_12_24) + (1|set) + log(n_genes))

saveRDS(res.proc, file = "braak_resproc_subtype_dx.RDS")

contrasts = c(Diff = 'braak_4_12_24AD - braak_4_12_24CTR')

res.dl = dreamlet( res.proc, ~ scale(age_4_12_24) + (1|gender_4_12_24) + (1|set) + log(n_genes) + braak_4_12_24 + 0, contrasts=contrasts)

saveRDS(res.dl, file = "braak_resdl_subtype_dx.RDS")

tab = topTable(res.dl, coef="Diff", number=Inf)

# compute FDR within each cell type
res = tab %>%
    as_tibble %>%
    group_by(assay) %>%
    mutate(FDR.within = p.adjust(P.Value, "fdr"))

# counts results
table(res$adj.P.Val < 0.05)
table(res$FDR.within < 0.05)


write.table(res, file="braak_subtype_dx.csv", sep=",", quote=FALSE)


go.gs = get_GeneOntology(to="SYMBOL")
res_zenith = zenith_gsa(res.dl, coef = 'Diff', go.gs)
write.table(res_zenith, file="braak_zenith_subtype_dx.csv", sep=",", row.names = FALSE, quote=FALSE)



res.proc = processAssays( pb,  ~ cdr_4_12_24 + scale(age_4_12_24) + (1|gender_4_12_24) + (1|set) + log(n_genes))

saveRDS(res.proc, file = "cdr_resproc_subtype_dx.RDS")

contrasts = c(Diff = 'cdr_4_12_24DEM - cdr_4_12_24CTR')

res.dl = dreamlet( res.proc, ~ scale(age_4_12_24) + (1|gender_4_12_24) + (1|set) + log(n_genes) + cdr_4_12_24 + 0, contrasts=contrasts)

saveRDS(res.dl, file = "cdr_resdl_subtype_dx.RDS")

tab = topTable(res.dl, coef="Diff", number=Inf)

# compute FDR within each cell type
res = tab %>%
    as_tibble %>%
    group_by(assay) %>%
    mutate(FDR.within = p.adjust(P.Value, "fdr"))

# counts results
table(res$adj.P.Val < 0.05)
table(res$FDR.within < 0.05)


write.table(res, file="cdr_subtype_dx.csv", sep=",", quote=FALSE)


go.gs = get_GeneOntology(to="SYMBOL")
res_zenith = zenith_gsa(res.dl, coef = 'Diff', go.gs)
write.table(res_zenith, file="cdr_zenith_subtype_dx.csv", sep=",", row.names = FALSE, quote=FALSE)

