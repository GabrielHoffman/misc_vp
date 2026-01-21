
# Differential expression based on somatic rate

# cd /hpc/users/hoffmg01/work/somatic_rate

suppressPackageStartupMessages({
library(tidyverse)
library(lme4)
library(variancePartition)
library(ggplot2)
library(dreamlet)
})

# Load estimated somatic rate
df = readRDS("PsychAD_SMrate_poolID.rds") %>%
      as_tibble

# Evaluate correlation across technical replicates
res = lapply( unique(df$CellType), function(CT){
  data = df %>%
    filter(CellType == CT)

  fit = lmer( SomaticRate ~ (1|SubID), data )  
  vp = calcVarPart(fit)

  data.frame(CellType = CT, R2_SubID = vp[[1]])
}) %>%
  bind_rows

res %>%
  ggplot(aes(R2_SubID, CellType)) +
    geom_bar(stat="identity") + 
    theme_classic() + 
    theme(aspect.ratio=1) +
    scale_x_continuous(lim=c(0,1), expand=c(0,0)) +
    xlab("Intra class correlation")




# Differential expression
##########################

# read Pseudobulk
file = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/pseudobulk/FULL_2024-02-01_18_49_PB_SubID_class.RDS"
pb = readRDS(file)

# merge with SomaticRate for each SubID and cell class
metadata(pb)$aggr_means = metadata(pb)$aggr_means %>%
  left_join(df %>% 
    dplyr::rename(class = CellType) %>%
    select(class, SubID, SomaticRate, Age) %>%
    group_by(class, SubID) %>%
    summarize(SomaticRate = mean(SomaticRate)) %>%
    distinct,
    by = c("class", "SubID"))

metadata(pb)$aggr_means$SomaticRate[1:3]

# Account for AD, SCZ, PD, Dementia
X = with(colData(pb), data.frame(c11x,c10x,c33x,r09x))

res = unlist(apply(X, 1, function(x){
  if( all(is.na(x)) ) return("NA")
  if( all( x[!is.na(x)] == "Control") ) return("Control")
  if( !is.na(x['c11x']) && x['c11x'] == "AD" ) return("AD")
  if( !is.na(x['c10x']) && x['c10x'] == "SCZ" ) return("SCZ")
  if( !is.na(x['c33x']) ) return(x['c33x'])
  if( !is.na(x['r09x'])) return(x['r09x'])
  return(NA)
  } ))

pb$DiseaseStatus = res


# voom-style normalization
formula = ~ scale(Age) + Sex + scale(PMI) + log(n_genes) + percent_mito + mito_genes + ribo_genes + mito_ribo + DiseaseStatus

res.proc <- processAssays(pb, formula)

formula = update(formula, ~ . + SomaticRate)

# run dreamlet
fit = dreamlet( res.proc, formula)

tab = topTable(fit, "SomaticRate", number=Inf) %>%
  as_tibble

tab %>%
  group_by(assay) %>%
  summarize(nFDR = sum(adj.P.Val < 0.05))

saveRDS( tab, file="somatic_rate_DE.RDS")



library(zenith)
go.gs <- get_GeneOntology("CC", to = "SYMBOL")

res_zenith <- zenith_gsa(fit, go.gs, "SomaticRate") %>%
  as_tibble

res_zenith %>%
  filter(NGenes < 1000) %>%
  arrange(FDR) %>%
  filter(assay == "Oligo") %>%
  data.frame %>%
  head(20)


saveRDS( res_zenith, file="somatic_rate_DE_genesets.RDS")


scp minerva:"/hpc/users/hoffmg01/work/somatic_rate/somat*" .








