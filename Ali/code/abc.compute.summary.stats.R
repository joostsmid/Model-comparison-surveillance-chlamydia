### R script to compute summary statistic-based discrepancy distance given a set of simulated population data from abc.simulate.chlamydia.R

## Summary statistic-based discrepancy distance function
compute.summary.stats <- function(simulated.population) {
  
  Nsim <- length(simulated.population[[1]]$output[,1,1])
  
  # Organised output from simulate.chlamydia into 'easy to understand' matrices of notifications and tests by age and sex
  mock.notifications.m <- array(0,c(Nsim,nyears,4))
  mock.notifications.m[,,1] <- simulated.population[[1]]$output[,35,8:(nyears+7)]
  mock.notifications.m[,,2] <- simulated.population[[2]]$output[,35,8:(nyears+7)]
  mock.notifications.m[,,3] <- simulated.population[[3]]$output[,35,8:(nyears+7)]
  mock.notifications.m[,,4] <- simulated.population[[4]]$output[,35,8:(nyears+7)]
  
  mock.notifications.f <- array(0,c(Nsim,nyears,4))
  mock.notifications.f[,,1] <- simulated.population[[5]]$output[,35,8:(nyears+7)]
  mock.notifications.f[,,2] <- simulated.population[[6]]$output[,35,8:(nyears+7)]
  mock.notifications.f[,,3] <- simulated.population[[7]]$output[,35,8:(nyears+7)]
  mock.notifications.f[,,4] <- simulated.population[[8]]$output[,35,8:(nyears+7)]

  
  mock.tested.m <- array(0,c(Nsim,nyears,4))
  mock.tested.m[,,1] <- simulated.population[[1]]$output[,36,8:(nyears+7)]
  mock.tested.m[,,2] <- simulated.population[[2]]$output[,36,8:(nyears+7)]
  mock.tested.m[,,3] <- simulated.population[[3]]$output[,36,8:(nyears+7)]
  mock.tested.m[,,4] <- simulated.population[[4]]$output[,36,8:(nyears+7)]
  
  mock.tested.f <- array(0,c(Nsim,nyears,4))
  mock.tested.f[,,1] <- simulated.population[[5]]$output[,36,8:(nyears+7)]
  mock.tested.f[,,2] <- simulated.population[[6]]$output[,36,8:(nyears+7)]
  mock.tested.f[,,3] <- simulated.population[[7]]$output[,36,8:(nyears+7)]
  mock.tested.f[,,4] <- simulated.population[[8]]$output[,36,8:(nyears+7)]
  
  # as discrepancy distance we take the sum of absolute fractional errors on each count (i.e. [ |mock - observed| / observed]) 
  # multiplied by a normalized weighting function that prioritizes the fit to our intermediate age cohorts (i.e. 15-19/20-24/25-29/30-34)
  # JS080618: the weighting function was deleted
  epsilon <- numeric(Nsim)
  for (i in 1:Nsim) {
    epsilon[i] <- sum(abs(notifications.m-t(mock.notifications.m[i,,]))/notifications.m) + 
    sum(abs((tested.m-t(mock.tested.m[i,,]))[!is.na(tested.m)])/tested.m[!is.na(tested.m)]) + 
    sum(abs(notifications.f-t(mock.notifications.f[i,,]))/notifications.f) + 
    sum(abs((tested.f-t(mock.tested.f[i,,]))[!is.na(tested.f)])/tested.f[!is.na(tested.f)]) 
    }
  
  return(epsilon)

}
