### R script to build mean vectors and covariance matrices with which to propose new parameters given our current ABC approximation to the posterior

build.mean.cov.matrices <- function(theta) {
  
  Nsim <- dim(theta)[1]
  
  ## Set up 'easy (or easier!) to understand' vectors of model infection and screening rates by age and sex, plus our intial set of 9 time-invariant parameters: we will then build (sample) mean and covariance matrices separately for each of these
  theta.fixed <- theta[,1:9]
  
  theta.m.inf <- theta[,(9+1):(9+nyears*4)]
  theta.m.inf <- log(theta.m.inf/(1-theta.m.inf))
  
  theta.f.inf <- theta[,((9+nyears*4)+1):(9+nyears*4*2)]
  theta.f.inf <- log(theta.f.inf/(1-theta.f.inf))
  
  theta.m.screening <- theta[,(9+nyears*4*2+1):(9+nyears*4*3)]
  theta.m.screening <- log(theta.m.screening/(1-theta.m.screening))
  
  theta.f.screening <- theta[,(9+nyears*4*3+1):(9+nyears*4*4)]
  theta.f.screening <- log(theta.f.screening/(1-theta.f.screening))
  
  mean.cov.master <- list()
  
  mean.cov.master$means <- list()
  mean.cov.master$covs <- list()
  
  mean.cov.master$means[[1]] <- colMeans(theta.fixed)
  mean.cov.master$means[[2]] <- colMeans(theta.m.inf)
  mean.cov.master$means[[3]] <- colMeans(theta.f.inf)
  mean.cov.master$means[[4]] <- colMeans(theta.m.screening)
  mean.cov.master$means[[5]] <- colMeans(theta.f.screening)
  
  mean.cov.master$covs[[1]] <- cov(theta.fixed)
  mean.cov.master$covs[[2]] <- cov(theta.m.inf)
  mean.cov.master$covs[[3]] <- cov(theta.f.inf)
  mean.cov.master$covs[[4]] <- cov(theta.m.screening)
  mean.cov.master$covs[[5]] <- cov(theta.f.screening)
  
  return(mean.cov.master)

}

