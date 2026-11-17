
# Differential analysis on SCENIC AUCs
# Nov 17, 2026

library(SingleCellExperiment)
library(fastglmm)
library(BatchRegression)
library(variancePartition)
library(parallel)
library(tidyverse)

file <- "/sc/arion/projects/CommonMind/seon/BD2/SCENIC_GRN/Combined_GRN_out/h5ad/combinedACC_AUCellLogitScores_noRawRNA_compressed_recode_gzip.h5ad"
sce <- GenomicDataStream::readH5AD(file, ondisk=FALSE)
# modify sce to account for Gaussian variable
cnt <- as.matrix(assay(sce, "counts"))
shift <- 1 - min(rowMins(cnt))
assay(sce, "counts") <- cnt + shift


# fixed effect formula
# drop set, it has K=34 levels and time is cubic in K
form <- ~ 0 + Dx + Sex + scale(PMI) + scale(n_genes) + percent_mito + mito_ribo #+ set 

# for each cell type
fitList = mclapply(unique(sce$psychad_class), function(CT){

  # extract cell type CT
  sce2 = sce[,sce$psychad_class==CT]

  # create fixed effect design matrix
  design <- model.matrix(form, colData(sce2))

  # create random effects design matrix
  Z <- preprocess_indicator(sce2$normalized_final_donor)

  # fit regression model
  fit <- glmmFitResponses(counts(sce2), design, Z, family=gaussian(), detail=3)

  fit
  }, mc.cores=12)
names(fitList) <- unique(sce$psychad_class)

res = lapply(names(fitList), function(CT){
  results(fitList[[CT]], coef="DxBD - DxControl") %>%
    tibble(CellType = CT, .)
}) %>%
  bind_rows %>%
  mutate(FDR = p.adjust(p.value, "BH")) %>%
  arrange(p.value)

out = "/sc/arion/projects/CommonMind/hoffman/grn/de_test.tsv"
write_tsv(res, out)


# dream: multiple random effects
#########

# formula
form <- ~ 0 + Dx + Sex + scale(PMI) + scale(n_genes) + percent_mito + mito_ribo + (1|set) + (1|normalized_final_donor)

# for each cell type
fitList = mclapply(unique(sce$psychad_class), function(CT){

  # extract cell type CT
  sce2 = sce[,sce$psychad_class==CT]

  L <- makeContrastsDream(form, colData(sce2), contrasts = c(Dx = "DxBD - DxControl"))

  fit <- dream(counts(sce2), form, colData(sce2), L=L )
  fit <- eBayes(fit)
  fit
  }, mc.cores=12)
names(fitList) <- unique(sce$psychad_class)


res2 = lapply(names(fitList), function(CT){
  topTable(fitList[[CT]], coef="Dx", number=Inf) %>%
    rownames_to_column("ID") %>%
    tibble(CellType = CT, .) %>%
    select(-adj.P.Val)
}) %>%
  bind_rows %>%
  mutate(FDR = p.adjust(P.Value, "BH")) %>%
  arrange(P.Value)

out = "/sc/arion/projects/CommonMind/hoffman/grn/de_test2.tsv"
write_tsv(res2, out)





