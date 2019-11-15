# Gabriel Hoffman
#
# November 15, 2019

# Code to combine QC metrics into a single table

# cd /sc/orga/projects/CommonMind/data/symlinks_for_mondale/RNAseq_Ajeet/62431/Processed/RAPiD


library(data.table)
statsList = list()

# featureCounts
file = "featureCounts/62431.exon.geneID.txt.summary"

res = read.table(file, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
res$Software = "featureCounts"
colnames(res)[1:2] = c("metric", "Value")
statsList[['featureCounts']] = res

# STAR
file = 'star/62431.Log.final.out'
cmd = paste("cat ", file, " | grep '|' | sed 's/\t//g'")
res = fread(cmd=cmd, sep='|')
res = res[-c(1:3),]
res$Software = 'STAR'
colnames(res)[1:2] = c("metric", "Value")
res$Value = as.numeric(gsub("%", "", res$Value))
res$metric = gsub("%", "Perc", res$metric)
statsList[['STAR']] = res

# PICARD - AlignmentSummaryMetrics
file = "qc_metrics/62431.AlignmentSummaryMetrics"
cmd = paste("grep -v '#'", file)
res = t(fread( cmd=cmd)[3,])
res = data.frame(metric = rownames(res), Value = res, stringsAsFactors=FALSE)
res = res[!is.na(res$Value),]
res$Software = "PICARD"
statsList[['AlignmentSummaryMetrics']] = res

# PICARD -DupMetrics
file = "qc_metrics/62431.DupMetrics"
cmd = paste("grep -v '#'", file)
res = suppressWarnings(t(fread( cmd=cmd, nrows=4)))
res = data.frame(metric = rownames(res), Value = res, stringsAsFactors=FALSE)
res = res[!is.na(res$Value),]
res$Software = "PICARD.DupMetrics"
statsList[['DupMetrics']] = res

# PICARD -GcBiasSummaryMetrics
file = "qc_metrics/62431.GcBiasSummaryMetrics"
cmd = paste("grep -v '#'", file)
res = t(fread( cmd=cmd, nrows=4))
res = data.frame(metric = rownames(res), Value = res, stringsAsFactors=FALSE)
res = res[!is.na(res$Value),]
res$Software = "PICARD.GcBiasSummaryMetrics"
statsList[['GcBiasSummaryMetrics']] = res



# PICARD -InsertSizeMetrics
file = "qc_metrics/62431.InsertSizeMetrics"
cmd = paste("grep -v '#'", file)
res = suppressWarnings(t(fread( cmd=cmd, nrows=4)))
res = data.frame(metric = rownames(res), Value = res, stringsAsFactors=FALSE)
res = res[!is.na(res$Value),]
res$Software = "PICARD.InsertSizeMetrics"
statsList[['InsertSizeMetrics']] = res



# PICARD -RNASeqMetrics
file = "qc_metrics/62431.RNASeqMetrics"
cmd = paste("grep -v '#'", file)
res = suppressWarnings(t(fread( cmd=cmd, nrows=4)))
res = data.frame(metric = rownames(res), Value = res, stringsAsFactors=FALSE)
res = res[!is.na(res$Value),]
res$Software = "PICARD.RNASeqMetrics"
statsList[['RNASeqMetrics']] = res


# combine
statsTable = do.call("rbind", statsList)








