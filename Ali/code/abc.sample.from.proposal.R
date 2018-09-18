### R script to define a function for simulating from the current proposal, defined by its set of mean and covariance matrices (mean.cov.master) computed by abc.build.mean.cov.matrices.R (used during the SMC-ABC fit procedure)

## Required packages
load.library(tmvtnorm) # (available from online R archive: CRAN)
load.library(fields) # (available from online R archive: CRAN)

## Posterior simulation function
simulate.from.proposal <- function(Nsim,mean.cov.master) {
  
  # Nsim = number of draws required
  
  theta <- matrix(0,nrow=Nsim,ncol=(nyears*4*2*2+9))
  # 4 is the number of fixed age cohorts, 2 the number of sexes (obviously)
  # and the last factor of 2 is because we have 2 MVN models: one for infections and one for screening
  # plus there are 9 (assumed) time-invariant model parameters (e.g. the true positive rate)
  
  # draw the 9 (assumed) time-invariant parameters
  theta[,1:9] <- rtmvnorm(Nsim,mean=mean.cov.master$means[[1]],sigma=mean.cov.master$covs[[1]],lower=rep(0,9),upper=rep(1,9))
  
  # draw the MVN infection model
  sgp <- rmvnorm(Nsim,mean=mean.cov.master$means[[2]],sigma=mean.cov.master$covs[[2]])
  tsgp <- exp(sgp)/(1+exp(sgp))
  theta[,(9+1):(9+nyears*4)] <- tsgp # infection matrix for males
  sgp <- rmvnorm(Nsim,mean=mean.cov.master$means[[3]],sigma=mean.cov.master$covs[[3]])
  tsgp <- exp(sgp)/(1+exp(sgp))
  theta[,((9+nyears*4)+1):(9+nyears*4*2)] <- tsgp # infection matrix for females
  
  # draw the MVN screening model
  sgp <- rmvnorm(Nsim,mean=mean.cov.master$means[[4]],sigma=mean.cov.master$covs[[4]])
  tsgp <- exp(sgp)/(1+exp(sgp))
  theta[,(9+nyears*4*2+1):(9+nyears*4*3)] <- tsgp # screening matrix for males
  
  sgp <- rmvnorm(Nsim,mean=mean.cov.master$means[[5]],sigma=mean.cov.master$covs[[5]])
  tsgp <- exp(sgp)/(1+exp(sgp))
  theta[,(9+nyears*4*3+1):(9+nyears*4*4)] <- tsgp # screening matrix for females
  
  return(theta)

}

