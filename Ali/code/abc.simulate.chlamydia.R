### R script to simulate notification, count, and prevalence estimates for a mock population given a matrix of model parameters

simulate.chlamydia <- function(theta) {

  Nsim <- dim(theta)[1]
  
  ## Set up 'easy (or easier!) to understand' vectors of model infection and screening rates by age and sex
  theta.m.inf <- theta[,(9+1):(9+nyears*4)]
  p.inf.m.15.19 <- theta.m.inf[,(nyears*0+1):(nyears*1)]
  p.inf.m.20.24 <- theta.m.inf[,(nyears*1+1):(nyears*2)]
  p.inf.m.25.34 <- theta.m.inf[,(nyears*2+1):(nyears*3)]
  p.inf.m.35.44 <- theta.m.inf[,(nyears*3+1):(nyears*4)]
  
  theta.f.inf <- theta[,((9+nyears*4)+1):(9+nyears*4*2)]
  p.inf.f.15.19 <- theta.f.inf[,(nyears*0+1):(nyears*1)]
  p.inf.f.20.24 <- theta.f.inf[,(nyears*1+1):(nyears*2)]
  p.inf.f.25.34 <- theta.f.inf[,(nyears*2+1):(nyears*3)]
  p.inf.f.35.44 <- theta.f.inf[,(nyears*3+1):(nyears*4)]
  
  theta.m.screening <- theta[,(9+nyears*4*2+1):(9+nyears*4*3)]
  p.screening.m.15.19 <- theta.m.screening[,1:nyears]
  p.screening.m.20.24 <- theta.m.screening[,(nyears+1):(nyears*2)]
  p.screening.m.25.34 <- theta.m.screening[,(nyears*2+1):(nyears*3)]
  p.screening.m.35.44 <- theta.m.screening[,(nyears*3+1):(nyears*4)]
  
  theta.f.screening <- theta[,(9+nyears*4*3+1):(9+nyears*4*4)]
  p.screening.f.15.19 <- theta.f.screening[,1:nyears]
  p.screening.f.20.24 <- theta.f.screening[,(nyears+1):(nyears*2)]
  p.screening.f.25.34 <- theta.f.screening[,(nyears*2+1):(nyears*3)]
  p.screening.f.35.44 <- theta.f.screening[,(nyears*3+1):(nyears*4)]
  
  ## Simulate mock data for each age-sex cohort
  simulated.population <- list()
  # Start simulation seven years prior to observed data to build up the steady-state infected population
  
  for (i in 1:2) {
    sex <- c("m","f")[i]
    for (j in 1:4) {
      age <- c("15.19","20.24","25.34","35.44")[j]
      m <- (i-1)*4+j
      simulated.population[[m]] <- list()
      simulated.population[[m]]$name <- paste(sex,".",age,sep="")
      simulated.population[[m]]$output <- array(0,c(Nsim,36,(nyears+7))) # the plus 7 is because we must run our model under the first years' parameters for at least 7 'dummy' years to build up a steady state population of infectees
      
      # in each year, for each age-sex cohort [nrow=Nsim,ncol=36] gives the matrix with columns:
      
      # 1) n.inf.tot
      # 2) n.inf.tot.new
      # 3) n.inf.tot.old
      # 4) n.inf.new.symp
      # 5) n.inf.new.asymp
      # 6) n.inf.old.symp
      # 7) n.inf.old.asymp
      # 8) n.inf.old.symp.inherited.from.last.year.same.age.cohort
      # 9) n.inf.old.asymp.inherited.from.last.year.same.age.cohort
      # 10) n.inf.old.symp.inherited.from.last.year.lower.age.cohort
      # 11) n.inf.old.asymp.inherited.from.last.year.lower.age.cohort
      # 12) n.inf.tot.symp
      # 13) n.inf.tot.asymp
      # 14) n.inf.symp.attend
      # 15) n.inf.symp.attend.tested
      # 16) n.inf.symp.attend.tested.rep
      # 17) n.inf.symp.attend.tested.true.pos
      # 18) n.inf.symp.attend.tested.true.pos.rep
      # 19) n.inf.asymp.screened
      # 20) n.inf.asymp.screened.rep
      # 21) n.inf.asymp.screened.true.pos
      # 22) n.inf.asymp.screened.true.pos.rep
      # 23) n.tot.uninf
      # 24) n.uninf.screened
      # 25) n.uninf.screened.rep
      # 26) n.uninf.screened.false.pos
      # 27) n.uninf.screened.false.pos.rep
      # 28) n.inf.symp.uncured.by.end.of.year
      # 29) n.inf.symp.uncured.by.end.of.year.graduating
      # 30) n.inf.symp.uncured.by.end.of.year.not.graduating
      # 31) n.inf.asymp.uncured.by.end.of.year
      # 32) n.inf.asymp.uncured.by.end.of.year.graduating
      # 33) n.inf.asymp.uncured.by.end.of.year.not.graduating
      # 34) prevalence.est
      # 35) notifications reported
      # 36) tests reported
      
      ### the following loop is not easy to understand, except by going line by line and drawing the update directions
      ### basically, it starts by doing the book-keeping for how many infected patients of each type (symptomatic and asymptomatic) remain in the age-sex cohort under study uncured since last year, 
      ### and how many such uncured patients have aged into this group from the cohort below
      ### it then simulates the number of new infections given the infection rate for this age-sex cohort and pushes them through our model of the reporting network to simulate mock numbers of notifications
      ### and test counts in the NNDSS
      ### finally it finishes with some book-keeping for uncured infections and decides how many will age out of this cohort by the end of the year
      for (k in 1:(nyears+7)) {
        
        if (k!=1 & m!=1 & m!=5) {
          simulated.population[[m]]$output[,10,k] <- simulated.population[[m-1]]$output[,29,k-1]
          simulated.population[[m]]$output[,11,k] <- simulated.population[[m-1]]$output[,32,k-1]
          simulated.population[[m]]$output[,8,k] <- simulated.population[[m]]$output[,30,k-1]
          simulated.population[[m]]$output[,9,k] <- simulated.population[[m]]$output[,33,k-1]
          simulated.population[[m]]$output[,6,k] <- simulated.population[[m]]$output[,8,k] + simulated.population[[m]]$output[,10,k]
          simulated.population[[m]]$output[,7,k] <- simulated.population[[m]]$output[,9,k] + simulated.population[[m]]$output[,11,k]
          simulated.population[[m]]$output[,3,k] <- simulated.population[[m]]$output[,6,k] + simulated.population[[m]]$output[,7,k]
        } else {
          simulated.population[[m]]$output[,3,k] <- 0
          simulated.population[[m]]$output[,6,k] <- 0
          simulated.population[[m]]$output[,7,k] <- 0
          simulated.population[[m]]$output[,8,k] <- 0
          simulated.population[[m]]$output[,9,k] <- 0
          simulated.population[[m]]$output[,10,k] <- 0
          simulated.population[[m]]$output[,11,k] <- 0}
        
        n.available <- rep(get(simulated.population[[m]]$name)[max(k-7,1)],Nsim) - simulated.population[[m]]$output[,3,k]
        n.available[n.available < 0] <- 1
        simulated.population[[m]]$output[,2,k] <- rbinom(Nsim,n.available,get(paste("p.inf.",simulated.population[[m]]$name,sep=""))[,max(k-7,1)])
        
        if (sex=="m") {
          simulated.population[[m]]$output[,5,k] <- rbinom(Nsim,simulated.population[[m]]$output[,2,k],theta[,6])
        } else {
          simulated.population[[m]]$output[,5,k] <- rbinom(Nsim,simulated.population[[m]]$output[,2,k],theta[,7])}
        
        simulated.population[[m]]$output[,4,k] <- simulated.population[[m]]$output[,2,k] - simulated.population[[m]]$output[,5,k]
        simulated.population[[m]]$output[,4,k][simulated.population[[m]]$output[,4,k] < 0] <- 1
        simulated.population[[m]]$output[,1,k] <- simulated.population[[m]]$output[,2,k] + simulated.population[[m]]$output[,3,k]
        simulated.population[[m]]$output[,23,k] <- get(simulated.population[[m]]$name)[max(k-7,1)] - simulated.population[[m]]$output[,1,k]
        simulated.population[[m]]$output[,23,k][simulated.population[[m]]$output[,23,k] < 0] <- 1
        simulated.population[[m]]$output[,12,k] <- simulated.population[[m]]$output[,4,k] + simulated.population[[m]]$output[,6,k]
        simulated.population[[m]]$output[,13,k] <- simulated.population[[m]]$output[,5,k] + simulated.population[[m]]$output[,7,k]
        
        simulated.population[[m]]$output[,14,k] <- rbinom(Nsim,simulated.population[[m]]$output[,12,k],theta[,1])
        simulated.population[[m]]$output[,15,k] <- rbinom(Nsim,simulated.population[[m]]$output[,14,k],theta[,2])
        simulated.population[[m]]$output[,17,k] <- rbinom(Nsim,simulated.population[[m]]$output[,15,k],theta[,3])
        simulated.population[[m]]$output[,18,k] <- rbinom(Nsim,simulated.population[[m]]$output[,17,k],theta[,5])
        simulated.population[[m]]$output[,16,k] <- simulated.population[[m]]$output[,18,k] + rbinom(Nsim,simulated.population[[m]]$output[,15,k] - simulated.population[[m]]$output[,17,k],theta[,5])
        
        simulated.population[[m]]$output[,19,k] <- rbinom(Nsim,simulated.population[[m]]$output[,13,k],get(paste("p.screening.",simulated.population[[m]]$name,sep=""))[,max(k-7,1)])
        simulated.population[[m]]$output[,21,k] <- rbinom(Nsim,simulated.population[[m]]$output[,19,k],theta[,3])
        simulated.population[[m]]$output[,22,k] <- rbinom(Nsim,simulated.population[[m]]$output[,21,k],theta[,5])
        simulated.population[[m]]$output[,20,k] <- simulated.population[[m]]$output[,22,k] + rbinom(Nsim,simulated.population[[m]]$output[,19,k] - simulated.population[[m]]$output[,21,k],theta[,5])
        
        simulated.population[[m]]$output[,24,k] <- rbinom(Nsim,simulated.population[[m]]$output[,23,k],get(paste("p.screening.",simulated.population[[m]]$name,sep=""))[,max(k-7,1)])
        simulated.population[[m]]$output[,26,k] <- rbinom(Nsim,simulated.population[[m]]$output[,24,k],theta[,4])
        simulated.population[[m]]$output[,27,k] <- rbinom(Nsim,simulated.population[[m]]$output[,26,k],theta[,5])
        simulated.population[[m]]$output[,25,k] <- simulated.population[[m]]$output[,27,k] + rbinom(Nsim,simulated.population[[m]]$output[,24,k] - simulated.population[[m]]$output[,26,k],theta[,5])
        
        simulated.population[[m]]$output[,35,k] <- simulated.population[[m]]$output[,18,k] + simulated.population[[m]]$output[,22,k] + simulated.population[[m]]$output[,27,k]
        simulated.population[[m]]$output[,36,k] <- simulated.population[[m]]$output[,16,k] + simulated.population[[m]]$output[,20,k] + simulated.population[[m]]$output[,25,k]
        
        simulated.population[[m]]$output[,28,k] <- simulated.population[[m]]$output[,12,k] - simulated.population[[m]]$output[,17,k]
        simulated.population[[m]]$output[,28,k] <- simulated.population[[m]]$output[,28,k] - rbinom(Nsim,simulated.population[[m]]$output[,28,k],theta[,8])
        simulated.population[[m]]$output[,28,k] <- simulated.population[[m]]$output[,28,k] - rbinom(Nsim,simulated.population[[m]]$output[,28,k],theta[,9])
        
        simulated.population[[m]]$output[,31,k] <- simulated.population[[m]]$output[,13,k] - simulated.population[[m]]$output[,21,k]
        simulated.population[[m]]$output[,31,k] <- simulated.population[[m]]$output[,31,k] - rbinom(Nsim,simulated.population[[m]]$output[,31,k],theta[,8])
        simulated.population[[m]]$output[,31,k] <- simulated.population[[m]]$output[,31,k] - rbinom(Nsim,simulated.population[[m]]$output[,31,k],theta[,9])
        
        
        if (m!=1 & m!=5 & m!=4 & m!=8) {
          # compute uncured infections that go over to next year, to other age cohort
          if (m  %in% c(1,2,5,6)){
            prob.inf.transcohort <- 1/5
          } else {
            prob.inf.transcohort <- 1/10
          }
          simulated.population[[m]]$output[,29,k] <- rbinom(Nsim,simulated.population[[m]]$output[,29,k],prob.inf.transcohort)
          simulated.population[[m]]$output[,30,k] <- simulated.population[[m]]$output[,28,k] - simulated.population[[m]]$output[,29,k]
          simulated.population[[m]]$output[,32,k] <- rbinom(Nsim,simulated.population[[m]]$output[,31,k],prob.inf.transcohort)
          simulated.population[[m]]$output[,33,k] <- simulated.population[[m]]$output[,31,k] - simulated.population[[m]]$output[,32,k]
        } else if (m==1 | m==5) {
          simulated.population[[m]]$output[,29,k] <- simulated.population[[m]]$output[,28,k]
          simulated.population[[m]]$output[,30,k] <- 0
          simulated.population[[m]]$output[,32,k] <- simulated.population[[m]]$output[,31,k]
          simulated.population[[m]]$output[,33,k] <- 0
        } else {
          simulated.population[[m]]$output[,29,k] <- 0
          simulated.population[[m]]$output[,30,k] <- simulated.population[[m]]$output[,28,k]
          simulated.population[[m]]$output[,32,k] <- 0
          simulated.population[[m]]$output[,33,k] <- simulated.population[[m]]$output[,31,k]
        }
        
        simulated.population[[m]]$output[,34,k] <- (simulated.population[[m]]$output[,3,k]/2+(simulated.population[[m]]$output[,28,k]+simulated.population[[m]]$output[,31,k])/2)/get(simulated.population[[m]]$name)[max(k-7,1)]
        
      }
    }
  }
  
  return(simulated.population)

}

