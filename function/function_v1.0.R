#####################################################################################
## stanjm.prep.3Y2T: prepare data for stanfit and the following AUC
## Output: 
##    A list storing data and info used in stan and AUC
#####################################################################################
#####################################################################################
stanjm.prep.3Y2T <- function(longFormula, survFormula, longda, survda, assoc = "value"){
  ############################################################################
  ######### standata for longitudinal submodels
  ############################################################################
  stan.data <- list()
  stan.data$qnodes <- 15
  wq <- gaussKronrod(stan.data$qnodes)$wq  ## add one 0 weight, a trick to calculate log(h_02)
  sq <- gaussKronrod(stan.data$qnodes)$sq
  
  K <- length(longFormula)
  terms1 <- attr(terms(longFormula[[1]]), "term.labels")
  terms2 <- attr(terms(longFormula[[2]]), "term.labels")
  terms3 <- attr(terms(longFormula[[3]]), "term.labels")
  formX1 <- reformulate(terms1[!grepl(" | ", terms1)])
  formX2 <- reformulate(terms2[!grepl(" | ", terms2)])
  formX3 <- reformulate(terms3[!grepl(" | ", terms3)])
  
  Y1_var <- all.vars(terms(longFormula[[1]]))[1]
  time_var <- strsplit(terms1[grepl(" | ", terms1)], " | ")[[1]][1]
  id_var <- strsplit(terms1[grepl(" | ", terms1)], " | ")[[1]][3]
  
  ## beta location for derivative of time
  stan.data$indFixed <- sapply(list(formX1, formX2, formX3), function(x) which(all.vars(x) == time_var)) + 1
  stan.data$indRandom <- cumsum(rep(2, K))
  
  stan.data$K <- K
  stan.data$y1 <- longda[, all.vars(terms(longFormula[[1]]))[1]]
  stan.data$y2 <- longda[, all.vars(terms(longFormula[[2]]))[1]]
  stan.data$y3 <- longda[, all.vars(terms(longFormula[[3]]))[1]]
  stan.data$y1_X <- model.matrix(formX1, longda)
  stan.data$y2_X <- model.matrix(formX2, longda)
  stan.data$y3_X <- model.matrix(formX3, longda)
  stan.data$y_N <- rep(nrow(longda), K)
  stan.data$y_K <- with(stan.data, c(ncol(y1_X), ncol(y2_X), ncol(y3_X))) 
  
  stan.data$b_N <- nrow(survda)
  stan.data$b_KM <- rep(2, K)           
  stan.data$b_K <- sum(stan.data$b_KM)  
  stan.data$y1_Z <- rbind(1, longda[[time_var]])  ## !!! only handle random intercept & slope !!! ##
  stan.data$y2_Z <- rbind(1, longda[[time_var]])
  stan.data$y3_Z <- rbind(1, longda[[time_var]])
  stan.data$y1_Z_id <- longda[[id_var]]
  stan.data$y2_Z_id <- longda[[id_var]]
  stan.data$y3_Z_id <- longda[[id_var]]
  
  ############################################################################
  ######### standata for main event submodel
  ############################################################################
  T_var <- all.vars(survFormula)[1]
  delta_var <- all.vars(survFormula)[2]
  W_var <- attr(terms(survFormula), "term.labels")
  if(!length(W_var)) W_var <- "1"
  formW <- reformulate(W_var)
  survmf <- model.matrix(formW, survda)
  
  stan.data$e_K <- ncol(survmf) - 1 
  stan.data$a_K <- K
  stan.data$Npat <- nrow(survda)
  stan.data$Nevents <- sum(survda[[delta_var]] == 1)
  stan.data$Npat_times_qnodes <- with(stan.data, Npat*qnodes)
  stan.data$nrow_e_Xq <- with(stan.data, Nevents + Npat_times_qnodes)
  stan.data$e_times <- c(survda[T_var][survda[delta_var] == 1],
                         survda[[T_var]] %*% t(1+sq)/2) 
  GK_index <- c(which(survda[delta_var] == 1), rep(1:nrow(survda), length(sq))) 
  stan.data$e_Xq <- survmf[GK_index, -1, drop = F]
  
  stan.data$basehaz_df <- 6
  library(splines)
  knots1 <- quantile(survda[T_var][survda[delta_var] == 1], c(1/3,2/3))
  Boundary.knots1 <- c(0, max(survda[T_var]))
  stan.data$basehaz_X <- with(stan.data, bs(e_times, df=basehaz_df, knots = knots1,
                                            Boundary.knots = Boundary.knots1, intercept = T))
  
  stan.data$qwts <- c(survda[[T_var]] %*% t(wq)/2)  ## Gauss-Kronrod approximation
  
  stan.data$nrow_y_Xq <- rep(stan.data$nrow_e_Xq, K)
  stan.data$y1_Xq <- model.matrix(formX1, survda)[GK_index, ] 
  stan.data$y1_Xq[, time_var] <- stan.data$e_times
  stan.data$y2_Xq <- model.matrix(formX2, survda)[GK_index, ]
  stan.data$y2_Xq[, time_var] <- stan.data$e_times
  stan.data$y3_Xq <- model.matrix(formX3, survda)[GK_index, ]
  stan.data$y3_Xq[, time_var] <- stan.data$e_times
  
  stan.data$y1_Zq <- rbind(1, stan.data$e_times)
  stan.data$y2_Zq <- rbind(1, stan.data$e_times)
  stan.data$y3_Zq <- rbind(1, stan.data$e_times)
  
  stan.data$y1_Zq_id <- survda[[id_var]][GK_index]
  stan.data$y2_Zq_id <- survda[[id_var]][GK_index]
  stan.data$y3_Zq_id <- survda[[id_var]][GK_index]
  
  ############################################################################
  ######### standata for dependent censoring submodel
  ############################################################################
  stan.data$e_K2 <- ncol(survmf) - 1
  stan.data$a_K2 <- K
  stan.data$Npat2 <- nrow(survda)
  stan.data$Nevents2 <- sum(survda[[delta_var]] == 2)
  stan.data$qnodes2 <- length(sq)
  stan.data$Npat_times_qnodes2 <- with(stan.data, Npat2*qnodes2)
  stan.data$nrow_e_Xq2 <- with(stan.data, Nevents2 + Npat_times_qnodes2)
  stan.data$e_times2 <- c(survda[T_var][survda[delta_var] == 2], survda[[T_var]] %*% t(1+sq)/2)
  GK_index2 <- c(which(survda[delta_var] == 2), rep(1:nrow(survda), length(sq))) ## T_i & T_i*(1+sq)/2 index
  stan.data$e_Xq2 <- survmf[GK_index2, -1, drop = F]  
  
  stan.data$basehaz_df2 <- 6
  knots2 <- quantile(survda[T_var][survda[delta_var] == 2], c(1/3,2/3))
  Boundary.knots2 <- c(0, max(survda[T_var]))
  stan.data$basehaz_X2 <- with(stan.data, bs(e_times2, df=basehaz_df2, knots = knots2,
                                             Boundary.knots = Boundary.knots2, intercept = T))
  
  stan.data$qwts2 <- c(survda[[T_var]] %*% t(wq)/2)
  
  stan.data$nrow_y_Xq2 <- rep(stan.data$nrow_e_Xq2, K)
  stan.data$y1_Xq2 <- model.matrix(formX1, survda)[GK_index2, ] 
  stan.data$y1_Xq2[, time_var] <- stan.data$e_times2
  stan.data$y2_Xq2 <- model.matrix(formX2, survda)[GK_index2, ]
  stan.data$y2_Xq2[, time_var] <- stan.data$e_times2
  stan.data$y3_Xq2 <- model.matrix(formX3, survda)[GK_index2, ]
  stan.data$y3_Xq2[, time_var] <- stan.data$e_times2
  
  
  stan.data$y1_Zq2 <- rbind(1, stan.data$e_times2)
  stan.data$y2_Zq2 <- rbind(1, stan.data$e_times2)
  stan.data$y3_Zq2 <- rbind(1, stan.data$e_times2)
  
  stan.data$y1_Zq_id2 <- survda[[id_var]][GK_index2]
  stan.data$y2_Zq_id2 <- survda[[id_var]][GK_index2]
  stan.data$y3_Zq_id2 <- survda[[id_var]][GK_index2]
  
  ############################################################################
  ######### stan.initial values
  ############################################################################
  sub.long1 <- lme4::lmer(longFormula[[1]], data = longda)
  sub.long2 <- lme4::lmer(longFormula[[2]], data = longda)
  sub.long3 <- lme4::glmer(longFormula[[3]], data = longda, family = binomial)
  sub.surv1 <- survival::coxph(as.formula(paste0("Surv(", T_var, ",", delta_var,"==1) ~", formW)[2]), data = survda)
  sub.surv2 <- survival::coxph(as.formula(paste0("Surv(", T_var, ",", delta_var,"==2) ~", formW)[2]), data = survda)
  
  stan.init <- list()
  stan.init$y1_beta <- sub.long1@beta
  stan.init$y2_beta <- sub.long2@beta
  stan.init$y3_beta <- sub.long3@beta
  stan.init$y1_aux <- summary(sub.long1)$sigma
  stan.init$y2_aux <- summary(sub.long2)$sigma
  stan.init$b_cholesky <- tryCatch(t(chol(as.matrix(bdiag(VarCorr(sub.long1)$id, VarCorr(sub.long2)$id, VarCorr(sub.long3)$id)))),
                                   error = function(e) {NULL})
  stan.init$b_mat <- cbind(ranef(sub.long1)$id, ranef(sub.long2)$id, ranef(sub.long3)$id)
  stan.init$e_beta <- as.array(coef(sub.surv1))
  stan.init$e_beta2 <- as.array(coef(sub.surv2))
  
  return(list(stan.data = stan.data, stan.init = stan.init,
              Y1_var=Y1_var, time_var=time_var, id_var=id_var, T_var=T_var, delta_var=delta_var,
              formX1=formX1, formX2=formX2,formX3=formX3,formW=formW,
              Boundary.knots1=Boundary.knots1,Boundary.knots2=Boundary.knots2,knots1=knots1,knots2=knots2,
              longFormula = longFormula, survFormula = survFormula, assoc = assoc))
}


#############################################################################################
## aucDep: Estimate the proposed net AUC in the presence of dependent censoring.
## Input: 
##    stanfit: fitted JM using stan (raw stan or rstanarm)
##    standat: stan preparation data which stores model info
##    newdata: longitudinal information for the predicted subject
##    T.start: numeric scalar denoting the time point up to which longitudinal information is to be used
##    T.end: numeric scalar denoting the time point for which a prediction of the survival status is of interest
##    Yvar: name of variable of interest
##    cutpoints: cutpoints to estimate sensitivity & specificity, can be a vector
##    simulate: TRUE for proposed estimator; FALSE for empirical Bayes estimator [See Rizopoulos(2011)]
##    M: number of Monte Carlo samples for random effects
##    n.pt: num of perturbation used to estimate variance caused by imputed event indicator
## Output: 
##    Estimation of auc, sen, spe, including their SE's
##
## Algorithms are extensions of those in Rizopoulos, D. (2011). The codes below for aucDep and related functions
## are adapted from the codes in R packages JM.
##############################################################################################
			       
aucDep <- function (stanfit, standat, newdata, T.start, T.end,
                    Yvar = NULL, cutpoints, simulate = TRUE, M = 200, n.pt = 200) {

  ############################### Extract info ###############################
  ############################################################################  
  if(class(stanfit) == "stanfit") stanfit <- stanfit@sim  ## extract @sim if using the raw stan
  Y_var <- standat$Y1_var                                 ## biomarker variable
  Y_var <- ifelse(is.null(Yvar), Y_var, Yvar)             ## whether to replace as Yvar
  time_var <- standat$time_var                            ## measurements time variable
  id_var <- standat$id_var
  T_var <- standat$T_var                                  ## observed event time variable 
  delta_var <- standat$delta_var
  
  ###################### Prepare data & set required data ####################
  ############################################################################ 
  newdata2 <- newdata[newdata[[T_var]] > T.start, ]        ## select alive at T.start
  newdata2 <- newdata2[newdata2[[time_var]] <= T.start, ]  ## to extract biomarker infor upto T.start
  
  last.rows <- tapply(1:nrow(newdata2), factor(newdata2[[id_var]]), tail, n = 1L) ## drop unused levels
  req.data <- newdata2[last.rows, ]  ## required data for AUC (to accomodate posterior_survfit)
  req.data$v12 <- 0                  ## event prob for censoring v1/v2
  
  uniq.id <- req.data[[id_var]]      ## uniq.id for the subjects used to calculate auc
  n.id <- length(uniq.id)
  scorei <- req.data[[Y_var]]        ## use last measurement befor T.start as score
  Ti <- req.data[[T_var]]            ## observed event time
  deltai <- req.data[[delta_var]]    ## event indicator
  
  ################# Calculate v1 & v2 for censored subjects ##################
  ############################################################################   
  cens.id <- uniq.id[deltai != 1 & Ti <= T.end]    ## censored b/w T.start and T.end
  score.ij <- (scorei %*% t(rep(1, n.id)) > t(scorei %*% t(rep(1, n.id)))) + 0.5*(scorei %*% t(rep(1, n.id)) == t(scorei %*% t(rep(1, n.id))))
  if(!simulate){
    v12 <- rep(0, n.id)
    v12[uniq.id %in% cens.id] <- survfitJM_modify(stanfit, standat, newdata[newdata[[id_var]] %in% cens.id, ], id_var = id_var, 
                                              delta_var = delta_var, last.time = Ti[uniq.id %in% cens.id], 
                                              time.to.pred = T.end, simulate = simulate, M = M)$out
    omegai <- (Ti <= T.end & deltai == 1) + (Ti <= T.end & deltai != 1) * v12
    omega.ij <- omegai %*% t(1 - omegai)
    auc.est <- sum(score.ij * omega.ij) / (sum(omega.ij)-sum(diag(omega.ij)))
    sen.est <- sapply(cutpoints, function(x) sum((scorei > x)*omegai) / sum(omegai))
    spe.est <- sapply(cutpoints, function(x) sum((scorei <= x)*(1 - omegai)) / sum(1 - omegai))
    
    list("auc" = auc.est, "sens" = sen.est, "spec" = spe.est)
  } else {
    v12.matrix <- matrix(0, nrow = n.id, ncol = M)   ## matrix to store the MCMC sample for censored id
    v12.matrix[uniq.id %in% cens.id, ] <- survfitJM_modify(stanfit, standat, newdata[newdata[[id_var]] %in% cens.id, ], id_var = id_var, 
                                                         delta_var = delta_var, last.time = Ti[uniq.id %in% cens.id], 
                                                         time.to.pred = T.end, simulate = simulate, M = M)$out
    auc.sample <- rep(NA, M)                 ## use v1,v2 sample to calculate auc
    sen.sample <- matrix(NA, M, length(cutpoints))
    spe.sample <- matrix(NA, M, length(cutpoints))
    for(ii in 1:M){
      omegai.s <- (Ti <= T.end & deltai == 1) + (Ti <= T.end & deltai != 1) * v12.matrix[ ,ii]
      omega.ij.s <- omegai.s %*% t(1 - omegai.s)
      auc.sample[ii] <- sum(score.ij * omega.ij.s) / (sum(omega.ij.s)-sum(diag(omega.ij.s)))
      sen.sample[ii, ] <- sapply(cutpoints, function(x) sum((scorei > x)*omegai.s) / sum(omegai.s))
      spe.sample[ii, ] <- sapply(cutpoints, function(x) sum((scorei <= x)*(1 - omegai.s)) / sum(1 - omegai.s))
    } 
    auc.est <- mean(auc.sample)               ## estimate of auc using MCMC sample
    sen.est <- apply(sen.sample, 2, mean)
    spe.est <- apply(spe.sample, 2, mean)
    
    auc.se.mc <- sd(auc.sample)               ## SE part caused by estimated thetas
    sen.se.mc <- apply(sen.sample, 2, sd)
    spe.se.mc <- apply(spe.sample, 2, sd)
    
    ############# purtabation part of auc.se assuming true thetas ##############
    ############################################################################ 
    if(n.pt <= 1){
      auc.se.pt <- NA
      sen.se.pt <- rep(NA, length(cutpoints))  ## sensitivity
      spe.se.pt <- rep(NA, length(cutpoints))  ## specificity
    }else{
      auc.pt <- matrix(NA, nrow = n.pt, ncol = M)                ## matrix to store auc.pt, for each pt, we have M's c-index
      sen.pt <- array(NA, dim = c(n.pt, M, length(cutpoints)))
      spe.pt <- array(NA, dim = c(n.pt, M, length(cutpoints)))
      for(i.pt in 1:n.pt){
        w.exp <- rexp(n.id)
        for(ii in 1:M){
          omegai.s <- (Ti <= T.end & deltai == 1) + (Ti <= T.end & deltai != 1) * v12.matrix[ ,ii]
          omega.ij.s <- omegai.s %*% t(1 - omegai.s)
          omega.ij.pt <- omega.ij.s * (w.exp %*% t(w.exp))
          auc.pt[i.pt, ii] <- sum(score.ij * omega.ij.pt) / (sum(omega.ij.pt)-sum(diag(omega.ij.pt)))
          
          omegai.pt <- omegai.s * w.exp
          sen.pt[i.pt, ii, ] <- sapply(cutpoints, function(x) sum((scorei > x)*omegai.pt) / sum(omegai.pt))
          spe.pt[i.pt, ii, ] <- sapply(cutpoints, function(x) sum((scorei <= x)*(w.exp - omegai.pt)) / sum(w.exp - omegai.pt))  ## !!!
        } 
      }
      auc.se.pt <- sd(rowMeans(auc.pt))                        ## perturbation part of auc.se
      sen.se.pt <- apply(apply(sen.pt, c(1,3), mean), 2, sd)
      spe.se.pt <- apply(apply(spe.pt, c(1,3), mean), 2, sd)   ## perturbation part of spe.se
    }
    auc.se <- sqrt(auc.se.mc^2 + auc.se.pt^2)  ## se has 2 parts: sandwich part & perturbation part
    sen.se <- sqrt(sen.se.mc^2 + sen.se.pt^2)
    spe.se <- sqrt(spe.se.mc^2 + spe.se.pt^2)
	  auc.res <- c(auc.est, auc.se, auc.se.mc, auc.se.pt)
    sens.res <- cbind(sen.est, sen.se, sen.se.mc, sen.se.pt)
    spec.res <- cbind(spe.est, spe.se, spe.se.mc, spe.se.pt)
    names(auc.res) <- colnames(sens.res) <- colnames(spec.res) <- c("Est", "SE", "SE.mc", "SE.pt")
    sens.res <- cbind("cutpoint" = cutpoints, sens.res)
    spec.res <- cbind("cutpoint" = cutpoints, spec.res)
    list("auc" = auc.res[1:2], "sens" = sens.res[, 1:3], "spec" = spec.res[, 1:3])	
  }
}


#############################################################################################
## survfitJM: Predict conditional probability of main event for ALL subjects at the given time
## **************** In fact, this is 1-SurvProb, NOT survival prob as the source code
## Input: 
##    stanfit: fitted JM using stan (raw stan or rstanarm)
##    standat: stan preparation data which stores model info
##    newdata: longitudinal information for the predicted subject (only one subject is allowed)
##    last.time: the last time when the subject is known to be alive; if NULL then the recent measurement time
##    time.to.pred: time to predict
##    simulate: TRUE for proposed estimator; FALSE for empirical Bayes estimator [See Rizopoulos(2011)]
##    M: number of Monte Carlo samples
##    scale: a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm
## Output: conditional probability of main event prob under dependent and independent censoring
##############################################################################################
survfitJM_modify <- function (stanfit, standat, newdata, id_var = "id", delta_var = "delta",
                          last.time = NULL, time.to.pred,
                          simulate = FALSE, M = 200, scale = 1.6) {

  ########################### Extract model info #############################
  ############################################################################
  time_var <- standat$time_var                  
  samples <- lapply(stanfit$samples, function(i) matrix(unlist(i), ncol = stanfit$iter, byrow = T))
  samples <- lapply(samples, function(i) i[, seq(stanfit$warmup + 1, stanfit$iter, stanfit$thin)])
  samples <- do.call(cbind, samples)
  samples <- samples[, sample(1:ncol(samples), M)]                          
  
  parameterization <- standat$assoc
  indFixed <- standat$stan.data$indFixed
  indRandom <- standat$stan.data$indRandom
  
  ########################### Extract longitudinal info ######################
  ############################################################################
  longFormula <- standat$longFormula
  survFormula <- standat$survFormula
  formX1 <- standat$formX1              
  formX2 <- standat$formX2
  formX3 <- standat$formX3
  formZ <- reformulate(time_var)       
  y1 <- newdata[, all.vars(terms(longFormula[[1]]))[1]]   
  y2 <- newdata[, all.vars(terms(longFormula[[2]]))[1]]   
  y3 <- newdata[, all.vars(terms(longFormula[[3]]))[1]]    
  X1 <- model.matrix(formX1, newdata)      
  X2 <- model.matrix(formX2, newdata)      
  X3 <- model.matrix(formX3, newdata)      
  Z <- model.matrix(formZ, newdata)        
  
  id <- as.numeric(unclass(newdata[[id_var]])) 
  id <- match(id, unique(id))                  
  
  ########################### Extract surv info ##############################
  ############################################################################
  data.id <- newdata[!duplicated(id), ]   
  delta <- data.id[[delta_var]]
  
  idT <- data.id[[id_var]]
  idT <- match(idT, unique(idT))
  
  formW <- standat$formW                  ## !!! must have one covariate for now, or have to revise stan !!!
  W <- model.matrix(formW, data.id)[, -1, drop=F]  
  
  ########################## Extract time.to.predict #########################
  ############################################################################
  obs.times <- split(newdata[[time_var]], id)            
  last.time <- if(is.null(last.time)) sapply(obs.times, max) else last.time
  time.to.pred <- if(length(time.to.pred) == 1) rep(time.to.pred, length(last.time))
  if(any(last.time > time.to.pred) | length(last.time) != length(idT) | length(time.to.pred) != length(idT))
    stop("time.to.pred should be larger than last.time")
  
  ########################### Extract parameter info #########################
  ############################################################################
  ncx1 <- ncol(X1)
  ncx2 <- ncol(X2)
  ncx3 <- ncol(X3)
  ncz <- ncol(Z)*standat$stan.data$K  ## !!!
  ncww <- ncol(W)

  varnames <- stanfit[["fnames_oi"]]
  betas1 <- samples[substr(varnames,1,7)=="y1_beta", , drop = F]
  betas2 <- samples[substr(varnames,1,7)=="y2_beta", , drop = F]
  betas3 <- samples[substr(varnames,1,7)=="y3_beta", , drop = F]
  sigma1 <- samples[substr(varnames,1,6)=="y1_aux", , drop = F]
  sigma2 <- samples[substr(varnames,1,6)=="y2_aux", , drop = F]
  D <- samples[substr(varnames,1,10)=="b_cholesky", , drop = F]
  gammas1 <- samples[substr(varnames,1,7)=="e_beta[", , drop = F]
  gammas2 <- samples[substr(varnames,1,7)=="e_beta2", , drop = F]
  alphas1 <- samples[substr(varnames,1,7)=="a_beta[", , drop = F]
  alphas2 <- samples[substr(varnames,1,7)=="a_beta2", , drop = F]
  gammas.bs1 <- samples[substr(varnames,1,6)=="e_aux[", , drop = F]
  gammas.bs2 <- samples[substr(varnames,1,6)=="e_aux2", , drop = F]
  
  Q <- standat$stan.data$qnodes       
  
  environment(log.posterior.b_modify) <- environment(S.b_modify) <- environment(ModelMats_modify) <- environment()
  ############################################################################
  ### construct model matrices to calculate the survival functions
  n.tp <- length(time.to.pred)                
  survMats <- survMats.last <- vector("list", n.tp)
  for(i in 1:n.tp){
    survMats[[i]] <- ModelMats_modify(time.to.pred[i], time_var, data.id[i, ], formX1, formX2, formX3, formZ)
    survMats.last[[i]] <- ModelMats_modify(last.time[i], time_var, data.id[i, ], formX1, formX2, formX3, formZ)
  }
  
  ### calculate the Empirical Bayes estimates and their (scaled) variance; 
  modes.b <- matrix(0, n.tp, ncz)       ## modes.b for each subject to be predicted
  Vars.b <- vector("list", n.tp)
  betas1.new <- rowMeans(betas1); betas2.new <- rowMeans(betas2); betas3.new <- rowMeans(betas3)
  sigma1.new <- rowMeans(sigma1); sigma2.new <- rowMeans(sigma2)
  D.chol <- matrix(rowMeans(D), standat$stan.data$b_K)           ## lower tri
  D.new <- D.chol %*% t(D.chol)
  gammas1.new <- rowMeans(gammas1); gammas2.new <- rowMeans(gammas2)
  alphas1.new <- rowMeans(alphas1); alphas2.new <- rowMeans(alphas2)
  gammas.bs1.new <- rowMeans(gammas.bs1); gammas.bs2.new <- rowMeans(gammas.bs2)
  for(i in 1:n.tp){
    ff <- function (b, y1, y2, y3, Mats, CR, ii) -log.posterior.b_modify(b, y1, y2, y3, Mats, CR, ii)  ## use info upto last.time
    CR <- if(delta[i] == 2) TRUE else if(delta[i] == 0) FALSE
    opt <- optim(rep(0, ncz), ff, y1 = y1[id==i], y2 = y2[id==i], y3 = y3[id==i], Mats = survMats.last[[i]], CR = CR, ii = i,
                 method = "BFGS", hessian = TRUE)
    modes.b[i, ] <- opt$par
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }
  
  if(!simulate) {
    v12 <- rep(NA, n.tp)
    for(i in 1:n.tp){
      S.last <- S.b_modify(last.time[i], modes.b[i, ], survMats.last[[i]], i)  ## use longitudinal upto last.time
      S.pred <- S.b_modify(time.to.pred[i], modes.b[i, ], survMats[[i]], i)
      v12[i] <- 1 - S.pred/S.last
    }
    return(list("out" = v12))
  } else {
    out <- matrix(NA, nrow = M, ncol = n.tp)     ## output to store v1, v2 (M*n.tp)
    success.rate <- matrix(FALSE, nrow = M, ncol = n.tp)
    b.old <- b.new <- modes.b
    if (n.tp == 1) dim(b.old) <- dim(b.new) <- c(1, ncz)  
    for (m in 1:M) {
      # Step 1: simulate new parameter values
      betas1.new <- betas1[,m]; betas2.new <- betas2[,m]; betas3.new <- betas3[,m]
      sigma1.new <- sigma1[,m]; sigma2.new <- sigma2[,m]
      D.chol <- matrix(D[,m], standat$stan.data$b_K)  ## lower tri
      D.new <- D.chol %*% t(D.chol)
      gammas1.new <- gammas1[,m]; gammas2.new <- gammas2[,m]
      alphas1.new <- alphas1[,m]; alphas2.new <- alphas2[,m]
      gammas.bs1.new <- gammas.bs1[,m]; gammas.bs2.new <- gammas.bs2[,m]
      
      v12 <- rep(NA, n.tp)
      for(i in 1:n.tp){
        # Step 2: simulate new random effects values
        proposed.b <- rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
        dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
        dmvt.proposed <- dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
        CR <- if(delta[i] == 2) TRUE else if(delta[i] == 0) FALSE
        a <- min(exp(log.posterior.b_modify(proposed.b, y1[id==i], y2[id==i], y3[id==i], survMats.last[[i]], CR, i) + dmvt.old -
                       log.posterior.b_modify(b.old[i, ], y1[id==i], y2[id==i], y3[id==i], survMats.last[[i]], CR, i) - dmvt.proposed), 1)
        ind <- (runif(1) <= a)
        success.rate[m, i] <- ind
        if (!is.na(ind) && ind) b.new[i, ] <- proposed.b
        
        # Step 3: compute Pr(T1* <= t| theta.new, b.new)
        S.last <- S.b_modify(last.time[i], b.new[i, ], survMats.last[[i]], i)
        S.pred <- S.b_modify(time.to.pred[i], b.new[i, ], survMats[[i]], i)
        v12[i] <- 1 - S.pred/S.last
      }
      b.old <- b.new
      out[m, ] <- v12
    }
    return(list("out" = t(out), "success.rate" = colMeans(success.rate)))
  }
}


#############################################################################################
## log.posterior.b: Calculate log of joint p(T, delta, Y, b)
##    require extra theta.hat: not an arg in this function (!!!use the .new pars)
## Input: 
##    b: random effects
##    y: observed longitudinal biomarkers 
##    Mats: output from ModelMats() to approximate integral in survival function
##    CR: TRUE if dependent censor observed; FALSE if independent censor observed
## Output: joint probability
##############################################################################################
log.posterior.b_modify <- function (b, y1, y2, y3, Mats, CR = T, ii) {

  id.i <- id %in% ii
  idT.i <- idT %in% ii
  X1.i <- X1[id.i, , drop = FALSE]  
  X2.i <- X2[id.i, , drop = FALSE]
  X3.i <- X3[id.i, , drop = FALSE]
  Z.i <- Z[id.i, , drop = FALSE] 
  W.i <- W[idT.i, , drop = FALSE] 
  
  mu.y1 <- c(X1.i %*% betas1.new) + rowSums(Z.i * rep(b[1:2], each = nrow(Z.i))) ## !!! to be more general
  mu.y2 <- c(X2.i %*% betas2.new) + rowSums(Z.i * rep(b[3:4], each = nrow(Z.i)))
  mu.y3 <- c(X3.i %*% betas3.new) + rowSums(Z.i * rep(b[5:6], each = nrow(Z.i)))
  logNorm1 <- dnorm(y1, mu.y1, sigma1.new, TRUE)
  logNorm2 <- dnorm(y2, mu.y2, sigma2.new, TRUE)
  logNorm3 <- dbinom(y3, 1, 1/(1+exp(-mu.y3)), TRUE)
  log.p.yb <- sum(logNorm1)+sum(logNorm2)+sum(logNorm3)   ## log{Pr(Y|b)}
  
  log.p.b <- dmvnorm(b, rep(0, ncz), D.new, TRUE)         ## log{Pr(b)}
  
  st <- Mats$st         ## use the whole history of Y; st = t(1+sq)/2
  wq <- Mats$wq         ## Gauss-Kronrod weight
  P <- Mats$P           ## P = time/2
  X1s <- Mats$X1s       ## X history
  X2s <- Mats$X2s
  X3s <- Mats$X3s
  Zs <- Mats$Zs

  if (parameterization %in% c("value")){   ## Ys = m_i(s), trajectory
    Y1s <- as.vector(X1s %*% betas1.new + rowSums(Zs * rep(b[1:2], each = nrow(Zs))))  ## !!!
    Y2s <- as.vector(X2s %*% betas2.new + rowSums(Zs * rep(b[3:4], each = nrow(Zs))))
    Y3s <- as.vector(X3s %*% betas3.new + rowSums(Zs * rep(b[5:6], each = nrow(Zs))))
  }                    

  gamma.t.w1 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas1.new), 0)  ## eta.tw = gamma %*% w
  gamma.t.w2 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas2.new), 0)  ## eta.tw = gamma %*% w
  
  ## baseline hazard: method == "spline-PH-GH": JM_R_book P53
  Bspline1 <- bs(st, df=standat$stan.data$basehaz_df, knots = standat$knots1,
                 Boundary.knots = standat$Boundary.knots1, intercept = T)
  Bspline2 <- bs(st, df=standat$stan.data$basehaz_df2, knots = standat$knots2,
                 Boundary.knots = standat$Boundary.knots2, intercept = T)
  logS.st1 <- exp(c(Bspline1 %*% gammas.bs1.new) + cbind(Y1s,Y2s,Y3s) %*% alphas1.new)    ## time term: h01(st) + alphas1*m_i(st)
  logS.st2 <- exp(c(Bspline2 %*% gammas.bs2.new) + cbind(Y1s,Y2s,Y3s) %*% alphas2.new)    ## h02(st) + alphas2*m_i(st)
  log.survival1 <- -exp(gamma.t.w1) * P * sum(wq * logS.st1)
  log.survival2 <- -exp(gamma.t.w2) * P * sum(wq * logS.st2)
  if (!CR) {
    return(log.survival1 + log.survival2 + log.p.yb + log.p.b)
  } else{
    log.haz2 <- log(logS.st2[length(logS.st2)]) + gamma.t.w2 ## log of dependent censor hazard
    return(log.haz2 + log.survival1 + log.survival2 + log.p.yb + log.p.b)
  }
}


#############################################################################################
## S.b: Calculate survival prob S(t|biomarker history, thetas) for the main event (e.g. death) 
##    require extra theta.hat: not an arg in this function (!!!use the .new pars)
## Input: 
##    t: time to predict
##    b: random effects
##    Mats: output from ModelMats() to approximate integral in survival function
## Output: survival probability
##############################################################################################
S.b_modify <- function (t, b, Mats, ii) { 
  idT.i <- idT %in% ii
  st <- Mats$st
  wq <- Mats$wq
  P <- Mats$P
  X1s <- Mats$X1s
  X2s <- Mats$X2s
  X3s <- Mats$X3s
  Zs <- Mats$Zs

  if (parameterization %in% c("value")){   ## Ys = m_i(s), trajectory
    Y1s <- as.vector(X1s %*% betas1.new + rowSums(Zs * rep(b[1:2], each = nrow(Zs))))
    Y2s <- as.vector(X2s %*% betas2.new + rowSums(Zs * rep(b[3:4], each = nrow(Zs))))
    Y3s <- as.vector(X3s %*% betas3.new + rowSums(Zs * rep(b[5:6], each = nrow(Zs))))
    Y1s.deriv <- 0
  }                    
  
  W.i <- W[idT.i, , drop = FALSE]
  gamma.t.w1 <- ifelse(!is.null(W.i), as.vector(W.i %*% gammas1.new), 0)  ## eta.tw = gamma %*% w
  
  ## baseline hazard: method == "spline-PH-GH": JM_R_book P53
  Bspline1 <- bs(st, df=standat$stan.data$basehaz_df, knots = standat$knots1,
                 Boundary.knots = standat$Boundary.knots1, intercept = T)
  logS.st1 <- exp(c(Bspline1 %*% gammas.bs1.new) + cbind(Y1s,Y2s,Y3s) %*% alphas1.new)    ## time term: h01(st) + alphas1*m_i(st)
  log.survival1 <- -exp(gamma.t.w1) * P * sum(wq * logS.st1)
  return(exp(log.survival1))
}


#############################################################################################
## ModelMats: Prepare for Gauss-Kronrod to approximate integral in survival function, 
##    but not engage survival function
## Input: 
##    time: the integral upper bound, i.e., use longitudinal history upto "time" in survival
##    time_var: time variable for the longitudinal measurements
##    data.id: survival data to build time-dependent X & Z history
##    formX: formula of X in longitudinal sub-model, to build X history
##    formZ: formula of Z in longitudinal sub-model, to build Z history
## Output: wq, sq, Xs & Zs matrix
##############################################################################################
ModelMats_modify <- function (time, time_var, data.id, formX1, formX2, formX3, formZ) {
  id.GK <- rep(1, each = standat$stan.data$qnodes + 1)  ## add 1 for log(h_02)
  
  wq <- c(gaussKronrod(standat$stan.data$qnodes)$wq, 0)  ## add one 0 weight, a trick to calculate log(h_02)
  sq <- gaussKronrod(standat$stan.data$qnodes)$sq
  P <- time / 2
  st <- c(P * (sq + 1), time)    ## st = t(1+sq)/2, add "time" to calculate log(h_02)
  
  data.id2 <- data.id[id.GK, ]   ## data.id: not an argument in this function
  data.id2[[time_var]] <- st
  
  out <- list(st = st, wq = wq, P = P)
  out$X1s <- model.matrix(formX1, data.id2)
  out$X2s <- model.matrix(formX2, data.id2)
  out$X3s <- model.matrix(formX3, data.id2)
  out$Zs <- model.matrix(formZ, data.id2)
  out
}


#############################################################################################
## Help functions from R packages JM.
#############################################################################################
gaussKronrod <- function (k = 15) {
  sq <- c(-0.949107912342758524526189684047851,
          -0.741531185599394439863864773280788,
          -0.405845151377397166906606412076961,
           0,
           0.405845151377397166906606412076961,
           0.741531185599394439863864773280788,
           0.949107912342758524526189684047851,
          -0.991455371120812639206854697526329,
          -0.864864423359769072789712788640926,
          -0.586087235467691130294144838258730,
          -0.207784955007898467600689403773245,
           0.207784955007898467600689403773245,
           0.586087235467691130294144838258730,
           0.864864423359769072789712788640926,
           0.991455371120812639206854697526329)
  wq15 <- c(0.063092092629978553290700663189204,
            0.140653259715525918745189590510238,
            0.190350578064785409913256402421014,
            0.209482141084727828012999174891714,
            0.190350578064785409913256402421014,
            0.140653259715525918745189590510238,
            0.063092092629978553290700663189204,
            0.022935322010529224963732008058970,
            0.104790010322250183839876322541518,
            0.169004726639267902826583426598550,
            0.204432940075298892414161999234649,
            0.204432940075298892414161999234649,
            0.169004726639267902826583426598550,
            0.104790010322250183839876322541518,
            0.022935322010529224963732008058970)
  wq7 <-  c(0.129484966168869693270611432679082,
            0.279705391489276667901467771423780,
            0.381830050505118944950369775488975,
            0.417959183673469387755102040816327,
            0.381830050505118944950369775488975,
            0.279705391489276667901467771423780,
            0.129484966168869693270611432679082)
  if (k == 7)
    list(sq = sq[1:7], wq = wq7)
  else
    list(sq = sq, wq = wq15)
}

dmvnorm <- function (x, mu, Sigma, log = FALSE) {
  if (!is.matrix(x))
    x <- rbind(x)
  p <- length(mu)
  if (p == 1) {
    dnorm(x, mu, sqrt(Sigma), log = log)
  } else {
    t1 <- length(mu) == length(Sigma)
    t2 <- all(abs(Sigma[lower.tri(Sigma)]) < sqrt(.Machine$double.eps))
    if (t1 || t2) {
      if (!t1)
        Sigma <- diag(Sigma)
      nx <- nrow(x)
      ff <- rowSums(dnorm(x, rep(mu, each = nx), sd = rep(sqrt(Sigma), each = nx), log = TRUE))
      if (log)
        ff
      else
        exp(ff)
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
      if (!all(ev >= -1e-06 * abs(ev[1])))
        stop("'Sigma' is not positive definite")
      ss <- x - rep(mu, each = nrow(x))
      inv.Sigma <- evec %*% (t(evec) / ev)
      quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
      fact <- -0.5 * (p * log(2 * pi) + sum(log(ev)))
      if (log)
        as.vector(fact - quad)
      else
        as.vector(exp(fact - quad))
    }
  }
}

rmvt <- function (n, mu, Sigma, df) {
  p <- length(mu)
  if (is.list(Sigma)) {
    ev <- Sigma$values
    evec <- Sigma$vectors
  } else {
    ed <- eigen(Sigma, symmetric = TRUE)
    ev <- ed$values
    evec <- ed$vectors
  }
  X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df) / df), each = p)
  if (n == 1L)
    drop(X)
  else
    t.default(X)
}

dmvt <- function (x, mu, Sigma, df, log = FALSE) {
  if (!is.numeric(x))
    stop("'x' must be a numeric matrix or vector")
  if (!is.matrix(x))
    x <- rbind(x)
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)) || ncol(x) != p)
    stop("incompatible arguments")
  ed <- eigen(Sigma, symmetric = TRUE)
  ev <- ed$values
  if (!all(ev >= -1e-06 * abs(ev[1])))
    stop("'Sigma' is not positive definite")
  ss <- x - rep(mu, each = nrow(x))
  inv.Sigma <- ed$vectors %*% (t(ed$vectors) / ev)
  quad <- rowSums((ss %*% inv.Sigma) * ss) / df
  fact <- lgamma((df + p) / 2) - lgamma(df / 2) - 0.5 * (p * (log(pi) + log(df)) + sum(log(ev)))
  if (log)
    fact - 0.5 * (df + p) * log(1 + quad)
  else
    exp(fact) * ((1 + quad) ^ (-(df + p) / 2))
}
