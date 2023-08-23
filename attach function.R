


# SIS for alpha
sis_alpha <- function(X, M, COV, p){
  s_alpha <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s1 <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s2 <- summary(fit)$coef[2]           #coefficients for alpha
    s3 <- summary(fit)$coef[2,4]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_alpha)))
  alpha_sis <- t(dat)
  colnames(alpha_sis) = colnames(M)
  return(s_alpha=alpha_sis)
}

# SIS for beta (the regression Y~X+M)
sis_beta1 <- function(X, M, Y, COV, p){
  s_beta <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(Y=Y, M = M[, j], X = X)
    } else {
      MX <- data.frame(Y=Y, M = M[, j], X = X, COV = COV)
    }
    fit <- coxph(Y ~., data = MX)
    s1 <- fit$var[1,1]                   #var for alpha
    s2 <- summary(fit)$coef[1]           #coefficients for alpha
    s3 <- summary(fit)$coef[1,5]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_beta)))
  beta_sis <- t(dat)
  colnames(beta_sis) = colnames(M)
  return(s_beta=beta_sis)
}

# SIS for beta (pcor for M and Y)
sis_beta2 <- function(X,M,Y,COV,p){
  p_cor = matrix(0,p,1)
  for(j in 1:p){
    if(is.null(COV)){
      d=data.frame(Y=Y[,1], M=M[,j], X=X)
      p_cor[j,] <- pcor(c(1,2,3),cov(d, method='spearman'))  #method of correlation
    }else{
      d=data.frame(Y=Y[,1], M=M[,j], X=X, COV=COV)            
      p_cor[j,] <- pcor(c(1,2,3,4,5), cov(d, method='spearman'))  #with 2 covariates
    }
  }
  rownames(p_cor) = colnames(M)
  return(p_cor = p_cor)
}

# Internal function: rdirichlet
# A function generate random number from Dirichlet distribution.

rdirichlet <- function (n = 1, alpha) 
{
  Gam <- matrix(0, n, length(alpha))
  for (i in 1:length(alpha)) Gam[, i] <- stats::rgamma(n, shape = alpha[i])
  Gam/rowSums(Gam)
}

JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}
nonnullPropEst <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) {
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}
nullParaEst<-function (x,gamma=0.1) 
{
  # x is a vector of z-values
  # gamma is a parameter, default is 0.1
  # output the estimated mean and standard deviation
  
  n = length(x)
  t = c(1:1000)/200
  
  gan    = n^(-gamma)
  that   = 0
  shat   = 0
  uhat   = 0
  epshat = 0
  
  phiplus   = rep(1,1000)
  phiminus  = rep(1,1000)
  dphiplus  = rep(1,1000)
  dphiminus = rep(1,1000)
  phi       = rep(1,1000)
  dphi      = rep(1,1000)
  
  for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind = min(c(1:1000)[(phi - gan) <= 0])
  tt = t[ind]
  a  = phiplus[ind]
  b  = phiminus[ind]
  da = dphiplus[ind]
  db = dphiminus[ind]
  c  = phi[ind]
  
  that   = tt
  shat   = -(a*da + b*db)/(tt*c*c)
  shat   = sqrt(shat)
  uhat   = -(da*b - db*a)/(c*c)
  epshat = 1 - c*exp((tt*shat)^2/2)
  
  return(musigma=list(mu=uhat,s=shat))
}

DACT = function(p_a,p_b,correction=NULL){
  Z_a = stats::qnorm(p_a,lower.tail = F)
  Z_b = stats::qnorm(p_b,lower.tail = F)
  pi0a = 1 - nonnullPropEst(Z_a,0,1)
  pi0b = 1 - nonnullPropEst(Z_b,0,1)
  # pi0a = locfdr::locfdr(Z_a,nulltype = 0)$fp0[5,3]
  # pi0b = locfdr::locfdr(Z_b,nulltype = 0)$fp0[5,3]
  if(pi0a > 1){
    pi0a = 1
  }
  if(pi0b >1){
    pi0b = 1
  }
  p.mat = cbind(p_a,p_b)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
  if(correction == "NULL"){
    p_dact = p_dact
  }
  if(correction == "Efron"){
    p_dact = EfronCorrect(p_dact)
  }
  if(correction == "JC"){
    p_dact = JCCorrect(p_dact)
  }
  return(p_dact)
}

