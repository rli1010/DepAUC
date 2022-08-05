## data split approach, with train size n=400 and test size n=400

## load required packages
packages <- c("lme4", "survival", "splines", "rstan", "JM")
suppressMessages(sapply(packages, require, character.only = TRUE))

## load required data and split as train and test set
load("Data.RData")
longtrain <- simda[simda$id <= 400, ]
longtest  <- simda[simda$id > 400, ]
survda    <- simda[tapply(1:nrow(simda), simda$id, tail, n = 1L), ]
survtrain <- survda[survda$id <= 400, ]

## prepare data for stan: Y1 ~ Normal, Y2 ~ Normal, Y3 ~ Bernoulli, current value association structure
source("function_v1.0.R")
stan.prep <- stanjm.prep.3Y2T(longFormula = list(Y1 ~ time+x1+(time|id), Y2 ~ time+x1+(time|id), Y3 ~ time+x1+(time|id)),
                              survFormula = Surv(obsT, delta) ~ x2, longda = longtrain, survda = survtrain, assoc = "value")

## fit joint model: took about 1 hr with iter = 2000; a fitted model was saved as "stanCRfit.RData"
stanCRfit <- stan(file = "jm_3Y2T.stan", data = stan.prep$stan.data, init = function() stan.prep$stan.init,
                  pars = "b_mat", include = FALSE, chains = 3, cores = 3, iter = 2000, seed = 123)
# load("stanCRfit.RData")
# traceplot(stanCRfit)

## auc estimation
## T.start, landmark time
## T.end, time point of interest for the main event
## Yvar, variable name of biomarker
## cutpoints, threshold in estimating sensitivity and specificity
## simulate = F to quickly obtain the point estimation, and also variance estimation if simulate = T
## M: number of MCMC resamples
## n.pt: number of perturbation in variance estimation
set.seed(123)
auc.net <- aucDep(stanCRfit, stan.prep, newdata = longtest, T.start = 3, T.end = 5, 
                  Yvar = "Y1", cutpoints = c(-0.5, -0.4), simulate = F, M = 200, n.pt = 200)
round(auc.net$auc, 3)
#   Est    SE 
# 0.613 0.040 





