### R script to define a function for returning the log density of the current proposal, defined by its set of mean and covariance matrices (mean.cov.master) computed by abc.build.mean.cov.matrices.R (used during the SMC-ABC fit procedure)

### Required packages
load.library(tmvtnorm) # (available from online R archive: CRAN)
load.library(fields) # (available from online R archive: CRAN)

## Log proposal density function
log.proposal.density <- function(theta,mean.cov.master) {
  
  Nsim <- dim(theta)[1]
  log.proposal.density <- numeric(Nsim)
  
  theta.fixed <- theta[,1:9]
  
  theta.m.inf <- theta[,(9+1):(9+nyears*4)]
  theta.m.inf <- log(theta.m.inf/(1-theta.m.inf))
  
  theta.f.inf <- theta[,((9+nyears*4)+1):(9+nyears*4*2)]
  theta.f.inf <- log(theta.f.inf/(1-theta.f.inf))
  
  theta.m.screening <- theta[,(9+nyears*4*2+1):(9+nyears*4*3)]
  theta.m.screening <- log(theta.m.screening/(1-theta.m.screening))
  
  theta.f.screening <- theta[,(9+nyears*4*3+1):(9+nyears*4*4)]
  theta.f.screening <- log(theta.f.screening/(1-theta.f.screening))
  
  log.proposal.density <- log.proposal.density + log(dtmvnorm(theta.fixed,mean=mean.cov.master$means[[1]],sigma=mean.cov.master$covs[[1]],lower=rep(0,9),upper=rep(1,9)))
  
  log.proposal.density <- log.proposal.density + dmvnorm(theta.m.inf,mean=mean.cov.master$means[[2]],sigma=mean.cov.master$covs[[2]],log=T)
  log.proposal.density <- log.proposal.density + dmvnorm(theta.f.inf,mean=mean.cov.master$means[[3]],sigma=mean.cov.master$covs[[3]],log=T)
  log.proposal.density <- log.proposal.density + dmvnorm(theta.m.screening,mean=mean.cov.master$means[[4]],sigma=mean.cov.master$covs[[4]],log=T)
  log.proposal.density <- log.proposal.density + dmvnorm(theta.f.screening,mean=mean.cov.master$means[[5]],sigma=mean.cov.master$covs[[5]],log=T)
  
  return(log.proposal.density)

}

