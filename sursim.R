
#' @examples
n=200            #sample size
p=1000         #dimension of mediators
alpha=rep(0,p)   #coefficients (mediator~exposure)
beta=rep(0,p)    #coefficients (outcome~mediators)
#' @param alpha a numeric vector specifying the regression coefficients alpha (exposure --> mediators).
#' @param beta a numeric vector specifying the regression coefficients beta (mediators --> outcome).
alpha[1:8] <- c(0.6,-0.5,0.65,0.45,0.5,0.4,0.45,0.7)
alpha[9] <- 0.5
beta[1:8] <- c(0.6,-0.5,0.45,0.6,0.65,0.7,0.55,0.3)
beta[10] <- 0.5
simsurBAMA = function(n, p, alpha, beta){
  set.seed(12345)
  X <- t(t(rbinom(n, 1, 0.6)))               #exposure
  gamma <-matrix(0.5,1,1)
  Z1 <- t(t(rbinom(n, 1, 0.3)))              #covariates Z1
  theta1 <- c(rep(0.3,p))                    #coefficients(Z1-->M)
  Z2 <- t(t(runif(n, 0, 1)))                 #covariates Z2
  theta2 <- c(rep(0.2,p))                    #coefficients(Z2-->M)
  Z <- cbind(Z1, Z2)
  phi <- c(0.3, -0.2)                        #coefficients(covariates-->outcome)
  ck <- t(runif(p, 0, 1))
  M <- matrix(0, n, p) 
  
  #mediators
  #i=1
  for(i in 1:n){
    e <- rnorm(p)
    M[i,] <- ck+X[i]*alpha+Z1[i]*theta1+Z2[i]*theta2+e
  }

  colnames(M) <- paste0("M", 1:ncol(M))
  haz <- 0.5*exp(0.5*X+0.3*Z[,1]-0.2*Z[,2]+0.6*M[,1]-0.5*M[,2]+0.45*M[,3]+0.6*M[,4]+0.65*M[,5]+
                   0.7*M[,6]+0.55*M[,7]+0.3*M[,8])   #baseline hazard function lambda0 <- 0.5
  ft <- rexp(n, haz)
  ct <- rexp(n, 2.5)              #censoring time
  time <- pmin(ft, ct)            #observed time
  status <- as.numeric(ft <= ct)  #censoring indicator
  Y <- survival::Surv(time, status)
  return(list(haz=haz, Y=Y, M=M, X=X, COV=Z, status=status,time=time,ft=ft))
  # return(list(haz=haz, Y=Y, M=M, X=X, COV=Z, status=status,time=time,ct=ct,ft=ft))
}

#' # Generate simulation data
#' simdat = sim_data(n, p, alpha, beta, seed=1029) 
sim_data<-simsurBAMA(n, p, alpha, beta)

