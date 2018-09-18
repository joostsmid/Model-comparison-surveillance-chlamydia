### R script to run the simulation algorithm a few more times to make posterior estimates of notification and test counts, incidence counts and fractions, positivity, & prevalence: 
### we could save these during runtime of the original SMC-ABC algorithm (abc.run.abc.R) but that's a bit of a waste of time in my opinion

# User specifications --------------------------------------------------------- 

# User specified utput folder and file for processing
fileNumber <- 65                     # User specified file number

# Run main script ---------------------------------------------------------------

# Load universal functions
source("code/load.library.R")

### load key modules
source("code/abc.read.in.data.R") # this will generate a report of 16 "errors" which may be ignored (these errors are due to the missing data in the 2007-2008 test counts)
source("code/abc.read.in.hyperparameters.R")
source("code/abc.simulate.chlamydia.R")
source("code/abc.compute.summary.stats.R")

# Load output file
outputFolder <- file.path(getwd(),"output")

load(file.path(outputFolder, paste("theta.test.",toString(fileNumber),".dat",sep = "")))

# Script parameters     
epsilon.thresh <- median(epsilon.current) 
Nsim <- length(epsilon.current)

simulated.population <- simulate.chlamydia(theta)
epsilon <- compute.summary.stats(simulated.population)

mock.not.m <- array(0,c(Nsim,nyears,4))
mock.not.m[,,1] <- simulated.population[[1]]$output[,35,8:(nyears+7)]
mock.not.m[,,2] <- simulated.population[[2]]$output[,35,8:(nyears+7)]
mock.not.m[,,3] <- simulated.population[[3]]$output[,35,8:(nyears+7)]
mock.not.m[,,4] <- simulated.population[[4]]$output[,35,8:(nyears+7)]

mock.not.f <- array(0,c(Nsim,nyears,4))
mock.not.f[,,1] <- simulated.population[[5]]$output[,35,8:(nyears+7)]
mock.not.f[,,2] <- simulated.population[[6]]$output[,35,8:(nyears+7)]
mock.not.f[,,3] <- simulated.population[[7]]$output[,35,8:(nyears+7)]
mock.not.f[,,4] <- simulated.population[[8]]$output[,35,8:(nyears+7)]

mock.test.m <- array(0,c(Nsim,nyears,4))
mock.test.m[,,1] <- simulated.population[[1]]$output[,36,8:(nyears+7)]
mock.test.m[,,2] <- simulated.population[[2]]$output[,36,8:(nyears+7)]
mock.test.m[,,3] <- simulated.population[[3]]$output[,36,8:(nyears+7)]
mock.test.m[,,4] <- simulated.population[[4]]$output[,36,8:(nyears+7)]

mock.test.f <- array(0,c(Nsim,nyears,4))
mock.test.f[,,1] <- simulated.population[[5]]$output[,36,8:(nyears+7)]
mock.test.f[,,2] <- simulated.population[[6]]$output[,36,8:(nyears+7)]
mock.test.f[,,3] <- simulated.population[[7]]$output[,36,8:(nyears+7)]
mock.test.f[,,4] <- simulated.population[[8]]$output[,36,8:(nyears+7)]

mock.inc.m <- array(0,c(Nsim,nyears,4))
mock.inc.m[,,1] <- simulated.population[[1]]$output[,2,8:(nyears+7)]
mock.inc.m[,,2] <- simulated.population[[2]]$output[,2,8:(nyears+7)]
mock.inc.m[,,3] <- simulated.population[[3]]$output[,2,8:(nyears+7)]
mock.inc.m[,,4] <- simulated.population[[4]]$output[,2,8:(nyears+7)]

mock.inc.f <- array(0,c(Nsim,nyears,4))
mock.inc.f[,,1] <- simulated.population[[5]]$output[,2,8:(nyears+7)]
mock.inc.f[,,2] <- simulated.population[[6]]$output[,2,8:(nyears+7)]
mock.inc.f[,,3] <- simulated.population[[7]]$output[,2,8:(nyears+7)]
mock.inc.f[,,4] <- simulated.population[[8]]$output[,2,8:(nyears+7)]

mock.incper.m <- array(0,c(Nsim,nyears,4))
mock.incper.m[,,1] <- simulated.population[[1]]$output[,2,8:(nyears+7)]/m.15.19
mock.incper.m[,,2] <- simulated.population[[2]]$output[,2,8:(nyears+7)]/m.20.24
mock.incper.m[,,3] <- simulated.population[[3]]$output[,2,8:(nyears+7)]/m.25.34
mock.incper.m[,,4] <- simulated.population[[4]]$output[,2,8:(nyears+7)]/m.35.44

mock.incper.f <- array(0,c(Nsim,nyears,4))
mock.incper.m[,,1] <- simulated.population[[5]]$output[,2,8:(nyears+7)]/f.15.19
mock.incper.m[,,2] <- simulated.population[[6]]$output[,2,8:(nyears+7)]/f.20.24
mock.incper.m[,,3] <- simulated.population[[7]]$output[,2,8:(nyears+7)]/f.25.34
mock.incper.m[,,4] <- simulated.population[[8]]$output[,2,8:(nyears+7)]/f.35.44

mock.prev.m <- array(0,c(Nsim,nyears,4))
mock.prev.m[,,1] <- simulated.population[[1]]$output[,34,8:(nyears+7)]
mock.prev.m[,,2] <- simulated.population[[2]]$output[,34,8:(nyears+7)]
mock.prev.m[,,3] <- simulated.population[[3]]$output[,34,8:(nyears+7)]
mock.prev.m[,,4] <- simulated.population[[4]]$output[,34,8:(nyears+7)]

mock.prev.f <- array(0,c(Nsim,nyears,4))
mock.prev.f[,,1] <- simulated.population[[5]]$output[,34,8:(nyears+7)]
mock.prev.f[,,2] <- simulated.population[[6]]$output[,34,8:(nyears+7)]
mock.prev.f[,,3] <- simulated.population[[7]]$output[,34,8:(nyears+7)]
mock.prev.f[,,4] <- simulated.population[[8]]$output[,34,8:(nyears+7)]

mock.pos.m <- array(0,c(Nsim,nyears,4))
mock.pos.m[,,1] <- simulated.population[[1]]$output[,35,8:(nyears+7)]/simulated.population[[1]]$output[,36,8:(nyears+7)]
mock.pos.m[,,2] <- simulated.population[[2]]$output[,35,8:(nyears+7)]/simulated.population[[2]]$output[,36,8:(nyears+7)]
mock.pos.m[,,3] <- simulated.population[[3]]$output[,35,8:(nyears+7)]/simulated.population[[3]]$output[,36,8:(nyears+7)]
mock.pos.m[,,4] <- simulated.population[[4]]$output[,35,8:(nyears+7)]/simulated.population[[4]]$output[,36,8:(nyears+7)]

mock.pos.f <- array(0,c(Nsim,nyears,4))
mock.pos.f[,,1] <- simulated.population[[5]]$output[,35,8:(nyears+7)]/simulated.population[[5]]$output[,36,8:(nyears+7)]
mock.pos.f[,,2] <- simulated.population[[6]]$output[,35,8:(nyears+7)]/simulated.population[[6]]$output[,36,8:(nyears+7)]
mock.pos.f[,,3] <- simulated.population[[7]]$output[,35,8:(nyears+7)]/simulated.population[[7]]$output[,36,8:(nyears+7)]
mock.pos.f[,,4] <- simulated.population[[8]]$output[,35,8:(nyears+7)]/simulated.population[[8]]$output[,36,8:(nyears+7)]

required <- numeric(Nsim)+1
required[epsilon < epsilon.thresh] <- 0

while (sum(required) > 6) { ### this loop may take a while ... it'll do the easy ones first but can take many tries to get the last few ... i sacrifice accuracy by allowing it to stop at 6 ... but if it's not being nice one might need to stop at 10 say (supposing Nsim ~ 10,000)
  
  simulated.population <- simulate.chlamydia(theta[which(required==1),])
  epsilon <- compute.summary.stats(simulated.population)
  
  if (length(which(epsilon < epsilon.thresh)) > 0) {
    
    valid.i <- which(epsilon < epsilon.thresh)
    valid.j <- which(required==1)[valid.i]
    
    mock.not.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,35,8:(nyears+7)]
    mock.not.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,35,8:(nyears+7)]
    mock.not.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,35,8:(nyears+7)]
    mock.not.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,35,8:(nyears+7)]
    
    mock.not.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,35,8:(nyears+7)]
    mock.not.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,35,8:(nyears+7)]
    mock.not.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,35,8:(nyears+7)]
    mock.not.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,35,8:(nyears+7)]
    
    
    mock.test.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,36,8:(nyears+7)]
    mock.test.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,36,8:(nyears+7)]
    mock.test.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,36,8:(nyears+7)]
    mock.test.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,36,8:(nyears+7)]
    
    mock.test.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,36,8:(nyears+7)]
    mock.test.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,36,8:(nyears+7)]
    mock.test.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,36,8:(nyears+7)]
    mock.test.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,36,8:(nyears+7)]
    
    mock.inc.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,2,8:(nyears+7)]
    
    mock.inc.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,2,8:(nyears+7)]
    mock.inc.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,2,8:(nyears+7)]
    

    mock.incper.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,2,8:(nyears+7)]/m.15.19
    mock.incper.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,2,8:(nyears+7)]/m.20.24
    mock.incper.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,2,8:(nyears+7)]/m.25.34
    mock.incper.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,2,8:(nyears+7)]/m.35.44
    
    mock.incper.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,2,8:(nyears+7)]/f.15.19
    mock.incper.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,2,8:(nyears+7)]/f.20.24
    mock.incper.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,2,8:(nyears+7)]/f.25.34
    mock.incper.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,2,8:(nyears+7)]/f.35.44
    
    mock.prev.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,34,8:(nyears+7)]
    
    mock.prev.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,34,8:(nyears+7)]
    mock.prev.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,34,8:(nyears+7)]
    
    

    mock.pos.m[valid.j,,1] <- simulated.population[[1]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.m[valid.j,,2] <- simulated.population[[2]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.m[valid.j,,3] <- simulated.population[[3]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.m[valid.j,,4] <- simulated.population[[4]]$output[valid.i,35,8:(nyears+7)]
    
    mock.pos.f[valid.j,,1] <- simulated.population[[5]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.f[valid.j,,2] <- simulated.population[[6]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.f[valid.j,,3] <- simulated.population[[7]]$output[valid.i,35,8:(nyears+7)]
    mock.pos.f[valid.j,,4] <- simulated.population[[8]]$output[valid.i,35,8:(nyears+7)]
    
    required[valid.j] <- 0
  }
  cat(sum(required),"\n")
}

# save the simplified output in a file called 'posterior.dat'
save(mock.not.m,mock.not.f,mock.test.m,mock.test.f,mock.inc.m,mock.inc.f,
     mock.incper.m,mock.incper.f,mock.prev.m,mock.prev.f,mock.pos.m,mock.pos.f,
     theta, file=file.path(outputFolder, paste0("posterior_",dataset,".dat"))) # "incidence.csv"