### R script to define a function for computing the logarithm of the prior density for a given matrix of model parameters

### Required packages
load.library(mvtnorm) # (available from online R archive: CRAN)
load.library(fields) # (available from online R archive: CRAN)

## Log prior density function
log.prior.density <- function(theta) {
  
  Nsim <- dim(theta)[1]
  log.prior.density <- numeric(Nsim)
  
  # log prior density contributions of the 9 (assumed) time-invariant parameters
  log.prior.density <- log.prior.density + log(dbeta(theta[,1],207,22))
  log.prior.density <- log.prior.density + log(dbeta(theta[,2],83,4))
  log.prior.density <- log.prior.density + log(dbeta(theta[,3],110,7))
  log.prior.density <- log.prior.density + log(dbeta(theta[,4],2,250))
  log.prior.density <- log.prior.density + log(dbeta(theta[,5],150,1))
  log.prior.density <- log.prior.density + log(dbeta(theta[,6],400,40))
  log.prior.density <- log.prior.density + log(dbeta(theta[,7],400,75))
  log.prior.density <- log.prior.density + log(dbeta(theta[,8],22,26))
  log.prior.density <- log.prior.density + log(dbeta(theta[,9],29,290))
  
  # log prior density contributions of the MVN infection model
  x <- 1:nyears
  y <- c(1,10,20,30) ### ?????????????????????????? y <- c(1,30,40,45,55,70) # the specific distance scaling here gives a correlation across age groups at fixed year that is significantly "weaker" than that at fixed age across the years
  z <- matrix(0,ncol=2,nrow=(4*nyears))
  for (i in 1:4) {z[((i-1)*nyears+1):(i*nyears),2] <- y[i]
                  z[((i-1)*nyears+1):(i*nyears),1] <- x}
  # apply the Matern correlation function to transform the distance matrix to a covariance matrix
  # p <- Matern(as.matrix(dist(z,diag=T,upper=T)),scale=prior.mvn.inf[2],range=prior.mvn.inf[3],smoothness=prior.mvn.inf[4])
  p <- prior.mvn.inf[2] * Matern(as.matrix(dist(z,diag=T,upper=T)),range=prior.mvn.inf[3],smoothness=prior.mvn.inf[4])
  theta.m.inf <- log(theta[,(9+1):(9+nyears*4)]/(1-theta[,(9+1):(9+nyears*4)]))
  log.prior.density <- log.prior.density + dmvnorm(theta.m.inf,mean=rep(prior.mvn.inf[1],nyears*4),sigma=p,log=T)
  theta.f.inf <- log(theta[,((9+nyears*4)+1):(9+nyears*4*2)]/(1-theta[,((9+nyears*4)+1):(9+nyears*4*2)]))
  log.prior.density <- log.prior.density + dmvnorm(theta.f.inf,mean=rep(prior.mvn.inf[1],nyears*4),sigma=p,log=T)
  
  # draw the MVN screening model
  x <- 1:nyears
  y <- c(1,10,20,30) ### ?????????????????????????? y <- c(1,30,40,45,55,70) # the specific distance scaling here gives a correlation across age groups at fixed year that is significantly "weaker" than that at fixed age across the years
  z <- matrix(0,ncol=2,nrow=(4*nyears))
  for (i in 1:4) {z[((i-1)*nyears+1):(i*nyears),2] <- y[i]
                  z[((i-1)*nyears+1):(i*nyears),1] <- x}
  # apply the Matern correlation function to transform the distance matrix to a covariance matrix
  # p <- Matern(as.matrix(dist(z,diag=T,upper=T)),scale=prior.mvn.screening[2],range=prior.mvn.screening[3],smoothness=prior.mvn.screening[4])
  p <- prior.mvn.screening[2] * Matern(as.matrix(dist(z,diag=T,upper=T)),range=prior.mvn.screening[3],smoothness=prior.mvn.screening[4])
  theta.m.screening <- log(theta[,(9+nyears*4*2+1):(9+nyears*4*3)]/(1-theta[,(9+nyears*4*2+1):(9+nyears*4*3)]))
  log.prior.density <- log.prior.density + dmvnorm(theta.m.screening,mean=rep(prior.mvn.screening[1],nyears*4),sigma=p,log=T)
  theta.f.screening <- log(theta[,(9+nyears*4*3+1):(9+nyears*4*4)]/(1-theta[,(9+nyears*4*3+1):(9+nyears*4*4)]))
  log.prior.density <- log.prior.density + dmvnorm(theta.f.screening,mean=rep(prior.mvn.screening[1],nyears*4),sigma=p,log=T)
  
  return(log.prior.density)

}

