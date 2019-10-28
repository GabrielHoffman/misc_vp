
## Variance partitioning for logistic regression

```r
library(variancePartition)

# load calcVarPart() for glm fit
# must load variancePartition first
source("https://raw.githubusercontent.com/GabrielHoffman/misc_vp/master/calcVarPart.R")

# geta dataset with two categories
data = iris[iris$Species %in% c("virginica", "versicolor"),]
data$Category = sample(c("A", "B", "C", "D"), nrow(data), replace=TRUE)

# fit logistic regression
form = Species ~  Petal.Length + Category 
fit.glm = glm(form, data, family=binomial())

# Run variance partitioning
calcVarPart(fit.glm)


# Fit Gaussian model forwith both functions
#-----------------------------------------

# Show values are the same

form = as.numeric(Species) ~ Petal.Length + Category
fit.lm = lm(form, data)

form = as.numeric(Species) ~  Petal.Length + Category
fit.glm = glm(form, data, family=gaussian())

calcVarPart(fit.lm)
calcVarPart(fit.glm)
```