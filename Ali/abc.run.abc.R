### R script to run the SMC-ABC algorithm for constraint of our model parameters

# Initilize

# choose dataset for UK
dataset <- "mindata"
#dataset <- "maxdata"

### load key modules
source("code/load.library.R")
source("code/abc.read.in.data.R") # this will generate a report of 16 "errors" which may be ignored (these errors are due to the missing data in the 2007-2008 test counts)
source("code/abc.read.in.hyperparameters.R")
source("code/abc.sample.from.prior.R")
source("code/abc.prior.log.density.R")
source("code/abc.simulate.chlamydia.R")
source("code/abc.compute.summary.stats.R")
source("code/abc.build.mean.cov.matrices.R")
source("code/abc.sample.from.proposal.R")
source("code/abc.proposal.log.density.R")

### Define 'output' folder
outputFolder <- file.path(getwd(),"output")

### SMC-ABC control parameters
Nsim <- 5000 # number of desired draws from the posterior
target.refreshment.rate <- 0.9 # target for number of unique particles in posterior approximation after running the MCMC refreshment kernal
discard.fraction <- 0.75 # what fraction of the current particle population to discard at the start of the round ... the closer to 1 the faster (but less stable) the convergence towards the true posterior
first.round.discard.fraction <- 0.975 # the fraction of particles to discard in the rejection ABC runs used to construct our very first posterior approximation
repeat.count <- 5 # number of times to run the MCMC refreshment kernel ... this must be user specified for first round, but is then predicted for subsequent rounds given the previous round's acceptance rate and our target refreshment rate
max.repeat.count <- 200 # maximum number of times to feasibly run the MCMC kernel ... this determines the stopping condition for our SMC-ABC analysis ... of course, the posterior approximation at each round is saved progressively so you can kill the procedure early if necessary without losing all your work

### simulate from prior for first round of ABC
for (i in 1:(1/(1-first.round.discard.fraction))) {
  theta <- simulate.from.prior(Nsim)
  simulated.population <- simulate.chlamydia(theta)
  epsilon.current <- compute.summary.stats(simulated.population)
  cat(i,"\n")
  if (i==1) {
    theta.save <- theta[sort.list(epsilon.current)[1:(Nsim*(1-first.round.discard.fraction))],]
    epsilon.current.save <- epsilon.current[sort.list(epsilon.current)[1:(Nsim*(1-first.round.discard.fraction))]]} else {
      theta.save <- rbind(theta.save,theta[sort.list(epsilon.current)[1:(Nsim*(1-first.round.discard.fraction))],])
      epsilon.current.save <- c(epsilon.current.save,epsilon.current[sort.list(epsilon.current)[1:(Nsim*(1-first.round.discard.fraction))]])
    }
}

theta <- theta.save
epsilon.current <- epsilon.current.save
output.thresh <- quantile(epsilon.current,1-discard.fraction)

count <- 1
## Run SMC-ABC until the efficiency of the MCMC refreshment kernel falls too low (i.e., when the expected number of repeat.count required to achieve the target.refreshment.rate exceeds max.repeat.count)
while (repeat.count < max.repeat.count) {
  cat("SMC ABC Round ",count,", Threshold = ",output.thresh,", N[kernel] = ",repeat.count,"\n") ## watch the threshold convergence towards zero (though it'll probably stop well before then) ... for the discrepancy distance currently defined the output.thresh here can be taken as a mean fractional predictive error when divided by the size of the parameter space (i.e., frac.err = output.thresh/(9+nyears*4*2+nyears*4*2)) ... should be able to get frac.err < 0.1 at least in a few hours of runtime (e.g., leave it overnight) on a standard laptop
  
  resampled <- sample(which(epsilon.current < output.thresh),Nsim,replace=T)
  theta <- theta[resampled,]
  log.theta.prior.density <- log.prior.density(theta)
  
  mean.cov <- build.mean.cov.matrices(theta)
  log.theta.proposal.density <- log.proposal.density(theta,mean.cov)
  
  n.accepted <- 0
  
  for (i in 1:repeat.count) {
    
    proposal.list <- simulate.from.proposal(Nsim,mean.cov)
    log.proposal.density.new <- log.proposal.density(proposal.list,mean.cov)
    log.proposal.prior.density <- log.prior.density(proposal.list)
    proposal.epsilon <- compute.summary.stats(simulate.chlamydia(proposal.list))
    
    accepted.moves <- which(proposal.epsilon < output.thresh & log(runif(Nsim,0,1)) < log.proposal.prior.density - log.theta.prior.density + log.theta.proposal.density - log.proposal.density.new)
    n.accepted <- n.accepted + length(accepted.moves)
    non.accepted.moves <- which(!(1:Nsim %in% accepted.moves))
    
    theta <- rbind(proposal.list[accepted.moves,],theta[non.accepted.moves,])
    log.theta.prior.density <- c(log.proposal.prior.density[accepted.moves],log.theta.prior.density[non.accepted.moves])
    epsilon.current <- c(proposal.epsilon[accepted.moves],epsilon.current[non.accepted.moves])}
  
  output.thresh <- quantile(epsilon.current,1-discard.fraction)
  p.acc <- n.accepted/repeat.count/Nsim
  repeat.count <- max(1,as.integer(log(1-target.refreshment.rate)/log(1-p.acc)))
  save(theta,epsilon.current,log.theta.prior.density,mean.cov,file=paste(outputFolder,"/theta.test.",count,".dat",sep=""))
  count <- count + 1
  
}

