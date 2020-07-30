

library(ggplot2)

colorCodes = c(
  GABA = '#66A61E',
  GLU = '#E6AB02',
  Olig = '#E7298A',
  MgAs = '#7570B3', 
  Other = "grey",
  GABAGLU = "dodgerblue")

df = read.table("adjust_data_add_sd", header=TRUE)
df$Cell_type = factor(df$Cell_type, names(colorCodes))
df$Regulatory_element = factor(df$Regulatory_element, names(colorCodes))

# Plot all results
ggplot(df, aes(eQTL_type, log2(OR), fill=Cell_type)) + geom_bar(stat='identity', position='dodge') + theme_bw(12) + theme(aspect.ratio=1) + facet_wrap(~Regulatory_element, nrow=1) + scale_fill_manual(values = colorCodes) + geom_errorbar( aes(ymin = log2(OR - sd), ymax = log2(OR + sd)), width=0,   position = position_dodge(width = 0.9))

# Plot results related to "Other"
df2 = df[df$Regulatory_element == "Other",]

ggplot(df2, aes(eQTL_type, -log2(OR), fill=Cell_type)) + geom_bar(stat='identity', position='dodge') + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = colorCodes) + geom_errorbar( aes(ymin = -log2(OR - sd), ymax = -log2(OR + sd)), width=0,   position = position_dodge(width = 0.9))


ggplot(df2, aes(eQTL_type, -log2(OR), fill=Cell_type)) + geom_bar(stat='identity', position='dodge') + theme_bw(12) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + scale_fill_manual(values = colorCodes) + geom_errorbar( aes(ymin = -log2(OR - sd), ymax = -log2(OR + sd)), width=0,   position = position_dodge(width = 0.9)) + facet_wrap(~Cell_type)





