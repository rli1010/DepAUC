

## load required packages
packages <- c("lme4", "survival", "splines", "rstan","JM")
suppressMessages(sapply(packages, require, character.only = TRUE))

load("Example.RData")

extraY_var <- "Y1"                    ## biomarker of interest
T.start = 0                           ## landmark time
T.end = 5                             ## time point of interest
M = 50                                 ## number of samples to use from the original MCMC samples, default = 200 for real data analysis
n.pt <- 50                             ## number of perturbation, default = 200

cutpoints <- median(longtest[[extraY_var]])  ## threshold to calculate sensitivity and specificity

source("function.R")
auc.net <- aucDep(stanCRfit, stan.prep, newdata=longtest, T.start, T.end, simulate=T, M=M, extraY_var = extraY_var, cutpoints=cutpoints, n.pt=n.pt)
print(auc.net$auc)






