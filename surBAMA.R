
#' Bayesian shrinkage estimation of high dimensional causal mediation effects in survival data
#' 

#' @param X a vector of exposure. 
#' @param Y a vector of outcome. Can be either continuous or binary (0-1).
#' @param M a data frame or matrix of high-dimensional mediators. Rows represent samples, columns 
#' represent variables.
#' @param COV a data frame or matrix of covariates dataset. Covariates specified here will not participate 
#' penalization. Default = \code{NULL}. 
#' #' @param OT a vector of observed failure times.
#' @param status a vector of censoring indicator (\code{status = 1}: uncensored; \code{status = 0}: censored)
#' @param k a integer of SIS method d = kn/log(n).Default = \code{1}.
#' @param FDRcut FDR cutoff applied to define and select significant mediators. Default = \code{0.05}. 
#' \function {BAMAss} is used to estimate and test high-dimensional mediation effects for survival data by spike and slab prior.
#' \function {BAMAhs} is used to estimate and test high-dimensional mediation effects for survival data by horseshoe prior.

library(survival)
library(glmnet)
library(dplyr)
library(brms)
library(bayestestR)
library(rstanarm)
library(ggm)

#######BAMA-ss (spike and slab prior)########
BAMAss <- function(X, COV, M, time, status,k=1,FDRcut = 0.05,verbose = FALSE){
  
  MZ <- cbind(M,COV,X)
  n <- length(X)
  p <- dim(M)[2]
  
  if(is.null(COV)){
    q <- 0
  }else{
    q <- dim(COV)[2]
  }
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Mediators screening ...", "     (", Sys.time(), ")")
  
  d_0 <- k*round(n/log(n))
  beta_SIS <- matrix(0,1,p)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZ_SIS <- MZ[,ID_S]
    fit <- survival::coxph(survival::Surv(time, status) ~ MZ_SIS)
    beta_SIS[i] <- fit$coefficients[1]
  }
  
  alpha_SIS <- matrix(0,1,p)
  alphavar_SIS <- matrix(0,1,p)
  alphase_SIS <- matrix(0,1,p)
  alphap_SIS <- matrix(0,1,p)
  XZ <- cbind(X,COV)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
    se_a <- ls.diag(fit_a)$std.err[2]
    var_a<-se_a^2
    p_a<-ls.print(fit_a)[["coef.table"]][[1]][2,4]
    alphavar_SIS[i] <- var_a
    alphap_SIS[i] <- p_a
  }
  
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[min(p, d_0)])
  
  d <- length(ID_SIS)
  if(verbose) message("        ", d, " mediators selected from the screening.")
  
  M_SIS <- M[, ID_SIS]
  
  message("Step 2: spike and slab estimates ...", "     (", Sys.time(), ")")
  
  ########### spike and slab estimation of beta ###########
  
  P_beta_SIS <- matrix(0,1,d)
  beta_est <- matrix(0,1,d)
  beta_SE <- matrix(0,1,d)
  MZ_SIS <- MZ[,c(ID_SIS, (p+1):(p+q+1))]
  MZ_SIS_1 <- t(t(MZ_SIS[,1]))
  # 
  
  ## spike and slab prior specification
  nvar=d+3
  s0=0.05
  ssde=set_prior("for (j in 1:nvar)
               target += log_sum_exp(log(1-gamma[j])+normal_lpdf(b[j]|0,s0*tau),
                                     log(gamma[j])+normal_lpdf(b[j]|0,tau))",check=F)+
    set_prior("target += inv_gamma_lpdf(tau | 1, 10)",check=F)+
    set_prior("for (j in 1:nvar)
               target += beta_lpdf(gamma[j]|1,1)",check=F)
  stanvars=stanvar(scode="real<lower=0> tau;",block="parameters")+
    stanvar(x=s0,name="s0",scode="real s0;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
    stanvar(scode="vector<lower=0,upper=1>[nvar] gamma;",block="parameters")
  
  if (is.null(COV)) {
    XM <- cbind(M_SIS, X,time,status)
    fitb=brm(formula = time | cens(1-status) ~ . ,
             data =XM, family = cox(link = "log"),
             prior=ssde,stanvars=stanvars,
             chain=3,iter=5000)
  } else { 
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(M_SIS, X,COV,time,status)
    fitb=brm(formula = time | cens(1-status) ~ . ,
             data =XM_COV, family = cox(link = "log"),
             prior=ssde,stanvars=stanvars,
             chain=3,iter=5000)
  }
  
  zz1<-summary(fitb)
  zz<-zz1[["fixed"]][c(1:d+1),]
  
  ## estimation of alpha
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  XZ <- cbind(X,COV)
  
  
  for (i in 1:d){
    fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a)/se_a
    P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
    alpha_SIS_est[i] <- est_a
    alpha_SIS_SE[i] <- se_a
  }
  
  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", names(ID_p_non), "\n")
  
  
  beta_SIS_est <- zz[,1]  # The non-zero ss estimators of beta
  beta_SIS_se<-zz[,2]
  var_beta <-(beta_SIS_se)^2 
  
  #########################################################################
  ################################ STEP 3 #################################
  #########################################################################
  message("Step 3: Multiple-testing procedure ...", "     (", Sys.time(), ")")
  # 
  
  ab_est <- alpha_SIS_est * beta_SIS_est   # the estimator of alpha*beta
  
  #######test
  #multiple-testing  procedure
  z.test<- abs((beta_SIS_est/zz[,2]))
  P_beta <- as.matrix(2*(1-pnorm(z.test,0,1)))
  P_beta <-ifelse(P_beta==0,1e-20,P_beta)
  
  
  P_alpha <-ifelse(P_alpha_SIS==0,1e-20,P_alpha_SIS)
  fdrcut <- DACT(t(P_alpha),P_beta,correction="NULL")
  ID_fdr <- which(abs(fdrcut) <=FDRcut)
  
  if (length(ID_fdr) > 0){
    alpha_est <- alpha_SIS_est[ID_fdr]
    alpha_se <- alpha_SIS_SE[ID_fdr]
    beta_est <- beta_SIS_est[ID_fdr]
    beta_se <- beta_SIS_se[ID_fdr]
    ab_est<-ab_est[ID_fdr]
    ID <- ID_SIS[ID_fdr]
    P_max <- abs(fdrcut[ID_fdr])
  }
  
  results<- data.frame(ID,alpha = alpha_est, alpha_se=alpha_se, beta = beta_est,beta_se=beta_se,
                       ab_est= ab_est, 
                       P=P_max,
                       check.names = FALSE)
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  return(list(results))
}

#######BAMA-hs (horseshoe prior )##########

BAMAhs <- function(X, COV, M, time, status,k=1,FDRcut = 0.05,verbose = FALSE){
  
  MZ <- cbind(M,COV,X)
  n <- length(X)
  p <- dim(M)[2]
  
  if(is.null(COV)){
    q <- 0
  }else{
    q <- dim(COV)[2]
  }
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Mediators screening ...", "     (", Sys.time(), ")")
  
  d_0 <- k*round(n/log(n))
  beta_SIS <- matrix(0,1,p)
  for (i in 1:p){
    ID_S <- c(i, (p+1):(p+q+1))
    MZ_SIS <- MZ[,ID_S]
    fit <- survival::coxph(survival::Surv(time, status) ~ MZ_SIS)
    beta_SIS[i] <- fit$coefficients[1]
  }
  
  alpha_SIS <- matrix(0,1,p)
  alphavar_SIS <- matrix(0,1,p)
  alphase_SIS <- matrix(0,1,p)
  alphap_SIS <- matrix(0,1,p)
  XZ <- cbind(X,COV)
  for (i in 1:p){
    fit_a  <- lsfit(XZ,M[,i],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
    se_a <- ls.diag(fit_a)$std.err[2]
    var_a<-se_a^2
    p_a<-ls.print(fit_a)[["coef.table"]][[1]][2,4]
    alphavar_SIS[i] <- var_a
    alphap_SIS[i] <- p_a
  }
  
  ab_SIS <- alpha_SIS*beta_SIS
  ID_SIS  <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[min(p, d_0)])
  
  d <- length(ID_SIS)
  if(verbose) message("        ", d, " mediators selected from the screening.")
  
  M_SIS <- M[, ID_SIS]
  
  # C_M <-colnames(XM)
  # d <- length(ID_SIS)
  
  message("Step 2: horse estimates ...", "     (", Sys.time(), ")")
  
  ## estimation of beta
  P_beta_SIS <- matrix(0,1,d)
  beta_est <- matrix(0,1,d)
  beta_SE <- matrix(0,1,d)
  MZ_SIS <- MZ[,c(ID_SIS, (p+1):(p+q+1))]
  MZ_SIS_1 <- t(t(MZ_SIS[,1]))
  
  prior<-set_prior(horseshoe(df = 3, scale_global = 1, df_global = 1,
                             scale_slab = 2,df_slab = 4, par_ratio = 0.001,
                             autoscale = TRUE))
  
  if (is.null(COV)) {
    XM <- cbind(M_SIS, X,time,status)
    fitb<-brm(formula = time | cens(1-status) ~ . ,
              data =XM, family =  cox(link = "log"),
              prior,
              warmup = 1000, iter = 5000, chains = 3,
              control = list(adapt_delta = 0.95))
  } else { 
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(M_SIS, X,COV,time,status)
    fitb<-brm(formula = time | cens(1-status) ~ . ,
              data =XM_COV, family = cox(link = "log"), 
              prior,
              warmup = 1000, iter = 5000, chains = 3,
              control = list(adapt_delta = 0.97))
  }
  
  zz1<-summary(fitb)
  zz<-zz1[["fixed"]][c(1:d+1),]
  
  ## estimation of alpha
  alpha_SIS_est <- matrix(0,1,d)
  alpha_SIS_SE <- matrix(0,1,d)
  P_alpha_SIS <- matrix(0,1,d)
  XZ <- cbind(X,COV)
  
  
  for (i in 1:d){
    fit_a  <- lsfit(XZ,M[,ID_SIS[i]],intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a)/se_a
    P_alpha_SIS[i] <- 2*(1-pnorm(sd_1,0,1))  ## the SIS for alpha
    alpha_SIS_est[i] <- est_a
    alpha_SIS_SE[i] <- se_a
  }
  
  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", names(ID_p_non), "\n")
  
  
  beta_SIS_est <- zz[,1]  # The non-zero ss estimators of beta
  beta_SIS_se<-zz[,2]
  var_beta <-(beta_SIS_se)^2 
  # ID_p <- ID_SIS[ID_p_non]  # The index of the ID of non-zero beta in the Cox regression
  
  #########################################################################
  ################################ STEP 3 #################################
  #########################################################################
  message("Step 3: Multiple-testing procedure ...", "     (", Sys.time(), ")")
  
  ab_est <- alpha_SIS_est * beta_SIS_est   # the estimator of alpha*beta
  
  #multiple-testing  procedure
  z.test<- abs((beta_SIS_est/zz[,2]))
  P_beta <- as.matrix(2*(1-pnorm(z.test,0,1)))
  P_beta <-ifelse(P_beta==0,1e-20,P_beta)
  
  P_alpha <-ifelse(P_alpha_SIS==0,1e-20,P_alpha_SIS)
  fdrcut <- DACT(t(P_alpha),P_beta,correction="NULL")
  ID_fdr <- which(fdrcut <=FDRcut)
  
  if (length(ID_fdr) > 0){
    alpha_est <- alpha_SIS_est[ID_fdr]
    alpha_se <- alpha_SIS_SE[ID_fdr]
    beta_est <- beta_SIS_est[ID_fdr]
    beta_se <- beta_SIS_se[ID_fdr]
    ab_est<-ab_est[ID_fdr]
    ID <- ID_SIS[ID_fdr]
    P_max <- abs(fdrcut[ID_fdr])
  }
  
  
  results<- data.frame(ID,alpha = alpha_est,alpha_se=alpha_se, beta = beta_est,beta_se=beta_se,
                       ab_est= ab_est, 
                       P=P_max,
                       check.names = FALSE)
  
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  return(list(results))
}
