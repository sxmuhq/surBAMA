
#' @examples
n=500            #sample size
p=10000         #dimension of mediators
alpha=rep(0,p)   #coefficients (mediator~exposure)
beta=rep(0,p)    #coefficients (outcome~mediators)
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
alpha[1:12] <- c(0.55,0.45,-0.4,-0.45,0.5,0.6,-0.4,-0.46,-0.40,0.5,0,0)
beta[1:12] <- c(0.52,-0.45,0.4,0.4,-0.54,0.6,-0.4,-0.5,0,0,0.4,-0.5)
simsurBAMA = function(n, p, alpha, beta,c0){
  set.seed(12345)
  X <- matrix(rbinom(n, 1, 0.6))             #exposure
  gamma <-matrix(0.5,1,1)
  Z1 <- matrix(rbinom(n, 1, 0.3))            #covariates Z1
  theta1 <- matrix((rep(0.3,p)))             #coefficients(Z1-->M)
  Z2 <- matrix(runif(n, 0, 1))               #covariates Z2
  theta2 <- matrix(rep(0.2,p))               #coefficients(Z2-->M)
  Z <- cbind(Z1, Z2)
  theta<-cbind(theta1,theta2)
  phi <- t(matrix(cbind(0.3,-0.2)))          #coefficients(covariates-->outcome)
    
  #mediators
  M <- matrix(0, n, p) 
  e <-matrix(rnorm(n*p),n)
  M <- X%*%(alpha) + Z%*%t(theta)+e
  colnames(M) <- paste0("M", 1:ncol(M))

  ##outcome
  MZX <- cbind(M,Z,X)
  beta_gamma <- cbind(beta,phi,gamma)
  haz <- 0.5*exp(MZX%*%t(beta_gamma))  #baseline hazard function lambda0 <- 0.5
  ft <- rexp(n, haz)
  ct <- rexp(n, c0)               #censoring time
  time <- pmin(ft, ct)            #observed time
  status <- as.numeric(ft <= ct)  #censoring indicator
  Y <- survival::Surv(time, status)
  return(list(haz=haz, Y=Y, M=M, X=X, COV=Z, status=status,time=time))
}

## Generate simulation data
sim_data<-simsurBAMA(n, p, alpha, beta,c0)

