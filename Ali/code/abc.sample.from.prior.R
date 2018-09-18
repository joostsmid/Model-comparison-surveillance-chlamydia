### R script to define a function for simulating from our prior

### Required packages
load.library(mvtnorm) # (available from online R archive: CRAN)
load.library(fields) # (available from online R archive: CRAN)

## Prior simulation function
simulate.from.prior <- function(Nsim) {

  # Nsim = number of draws required
  
  theta <- matrix(0,nrow=Nsim,ncol=(nyears*4*2*2+9))
  # 4 is the number of fixed age cohorts, 2 the number of sexes (obviously)
  # and the last factor of 2 is because we have 2 MVN models: one for infections and one for screening
  # plus there are 9 (assumed) time-invariant model parameters (e.g. the true positive rate)
  
  # draw the 9 (assumed) time-invariant parameters
  theta[,1] <- rbeta(Nsim,207,22) # attend|symp [0.863,0.905,0.939] # 4000 420 [0.896,0.905,0.913]
  theta[,2] <- rbeta(Nsim,83,4) # test|symp [0.901,0.957,0.987]
  theta[,3] <- rbeta(Nsim,110,7) # true.pos [0.891,0.943,0.975]
  theta[,4] <- rbeta(Nsim,2,250) # false.pos [0.001,0.007,0.022]
  theta[,5] <- rbeta(Nsim,150,1) # p.rep [0.976,0.995,0.999] # 1500 8 [0.990,0.995,0.998]
  theta[,6] <- rbeta(Nsim,400,40) # asymp.m [0.881,0.910,0.934]
  theta[,7] <- rbeta(Nsim,400,75) # asymp.f [0.808,0.843,0.873]
  theta[,8] <- rbeta(Nsim,22,26) # p.cured.after.year [0.321,0.458,0.599]
  theta[,9] <- rbeta(Nsim,29,290) # p.background.antibiotic [0.062,0.090,0.125] # 1050 10000 [0.090,0.095,0.101]
  
  # draw the MVN infection model
  x <- 1:nyears
  y <- c(1,10,20,30) ### ?????????????????????????? c(1,30,40,45,55,70) # the specific distance scaling here gives a correlation across age groups at fixed year that is significantly "weaker" than that at fixed age across the years
  z <- matrix(0,ncol=2,nrow=(4*nyears))
  for (i in 1:4) {
    z[((i-1)*nyears+1):(i*nyears),2] <- y[i]
    z[((i-1)*nyears+1):(i*nyears),1] <- x
  }
  
  # apply the Matern correlation function to transform the distance matrix to a covariance matrix. Note had to change how the Matern function works because of an update to 
  # the fields package
# p <- Matern(as.matrix(dist(z,diag=T,upper=T)),scale=prior.mvn.inf[2],range=prior.mvn.inf[3],smoothness=prior.mvn.inf[4])
  p <- prior.mvn.inf[2] * Matern(as.matrix(dist(z,diag=T,upper=T)), range=prior.mvn.inf[3],smoothness=prior.mvn.inf[4])
  sgp <- rmvnorm(Nsim,mean=rep(prior.mvn.inf[1],nyears*4),sigma=p)
  tsgp <- exp(sgp)/(1+exp(sgp)) # logistic transformation to a probability between 0 and 1
  theta[,(9+1):(9+nyears*4)] <- tsgp # infection matrix for males
  sgp <- rmvnorm(Nsim,mean=rep(prior.mvn.inf[1],nyears*4),sigma=p)
  tsgp <- exp(sgp)/(1+exp(sgp)) # logistic transformation to a probability between 0 and 1
  theta[,((9+nyears*4)+1):(9+nyears*4*2)] <- tsgp # infection matrix for females
  
  # draw the MVN screening model
  x <- 1:nyears
  y <- c(1,10,20,30) ### ?????????????????????????? y <- c(1,30,40,45,55,70) # the specific distance scaling here gives a correlation across age groups at fixed year that is significantly "weaker" than that at fixed age across the years
  z <- matrix(0,ncol=2,nrow=(4*nyears))
  for (i in 1:4) {
    z[((i-1)*nyears+1):(i*nyears),2] <- y[i]
    z[((i-1)*nyears+1):(i*nyears),1] <- x
  }
  
  # apply the Matern correlation function to transform the distance matrix to a covariance matrix. Note had to change how the Matern function works because of an update to 
  # the fields package
  # p <- Matern(as.matrix(dist(z,diag=T,upper=T)),scale=prior.mvn.screening[2],range=prior.mvn.screening[3],smoothness=prior.mvn.screening[4])
  p <- prior.mvn.screening[2] * Matern(as.matrix(dist(z,diag=T,upper=T)),range=prior.mvn.screening[3],smoothness=prior.mvn.screening[4])
  sgp <- rmvnorm(Nsim,mean=rep(prior.mvn.screening[1],nyears*4),sigma=p)
  tsgp <- exp(sgp)/(1+exp(sgp)) # logistic transformation to a probability between 0 and 1
  theta[,(9+nyears*4*2+1):(9+nyears*4*3)] <- tsgp # screening matrix for males
  sgp <- rmvnorm(Nsim,mean=rep(prior.mvn.screening[1],nyears*4),sigma=p)
  tsgp <- exp(sgp)/(1+exp(sgp)) # logistic transformation to a probability between 0 and 1
  theta[,(9+nyears*4*3+1):(9+nyears*4*4)] <- tsgp # screening matrix for females
  
  return(theta)

}

