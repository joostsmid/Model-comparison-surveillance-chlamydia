if(!require("rriskDistributions")) install.packages("rriskDistributions");library(rriskDistributions)
if(!require("deSolve")) install.packages("deSolve");library(deSolve)
if(!require("rootSolve")) install.packages("rootSolve");library( rootSolve)
if(!require("ggplot2")) install.packages("ggplot2");library(ggplot2)
if(!require("grid")) install.packages("grid");library(grid)

##############################
# read in data
##############################

population_f = read.csv('data/population_f.csv')
tests_f_max = read.csv('data/tests_f_max.csv')
tests_f_min = read.csv('data/tests_f_min.csv')
diagnoses_f_max = read.csv('data/diagnoses_f_max.csv')
diagnoses_f_min = read.csv('data/diagnoses_f_min.csv')

population_m = read.csv('data/population_m.csv')
tests_m_max = read.csv('data/tests_m_max.csv')
tests_m_min = read.csv('data/tests_m_min.csv')
diagnoses_m_max = read.csv('data/diagnoses_m_max.csv')
diagnoses_m_min = read.csv('data/diagnoses_m_min.csv')

#####################################################################
# sample from parameter distributions
#####################################################################

n_sample = 1000

# proportion sexually active men
p_active_m = array(NA, dim=c(1,2,n_sample))

# Beta parameters from analysis in R
p_active_m[1,1,] = rbeta(n_sample, 313.4784, 146.67142) # 16-19
p_active_m[1,2,] = rbeta(n_sample, 694.6619, 37.21263) # 20-24

# population sexually active in separate (16) years

pop_active_m = array(NA, dim=c(16,2,n_sample))

for (i in 1:16){
  for (j in 1:2){
    pop_active_m[i,j,] <- rbinom(n_sample, population_m[i,j+1], p_active_m[1,j,])
  }
}

# proportion sexually active women
p_active_f = array(NA, dim=c(1,2,n_sample))

# Beta parameters from analysis in R
p_active_f[1,1,] = rbeta(n_sample, 313.4784, 146.67142) # 16-19
p_active_f[1,2,] = rbeta(n_sample, 694.6619, 37.21263) # 20-24

# population sexually active in separate (16) years

pop_active_f = array(NA, dim=c(16,2,n_sample))

for (i in 1:16){
  for (j in 1:2){
    pop_active_f[i,j,] <- rbinom(n_sample, population_f[i,j+1], p_active_f[1,j,])
  }
}


# testing and diagnosis rates, per person per year
test_rate_f_max = array(NA, dim=c(16,2,n_sample))
diag_rate_f_max = array(NA, dim=c(16,2,n_sample))
test_rate_f_min = array(NA, dim=c(16,2,n_sample))
diag_rate_f_min = array(NA, dim=c(16,2,n_sample))
test_rate_m_max = array(NA, dim=c(16,2,n_sample))
diag_rate_m_max = array(NA, dim=c(16,2,n_sample))
test_rate_m_min = array(NA, dim=c(16,2,n_sample))
diag_rate_m_min = array(NA, dim=c(16,2,n_sample))

for (i in 1:16){
  for (j in 1:2){
    test_rate_m_max[i,j,] = rgamma(n_sample, tests_m_max[i,j+1]*population_m[i,j+1]/100, 1)/pop_active_m[i,j,]
    diag_rate_m_max[i,j,] = rgamma(n_sample, diagnoses_m_max[i,j+1]*population_m[i,j+1]/100000, 1)/pop_active_m[i,j,]
    test_rate_m_min[i,j,] = rgamma(n_sample, tests_m_min[i,j+1]*population_m[i,j+1]/100, 1)/pop_active_m[i,j,]
    diag_rate_m_min[i,j,] = rgamma(n_sample, diagnoses_m_min[i,j+1]*population_m[i,j+1]/100000, 1)/pop_active_m[i,j,]
    
    test_rate_f_max[i,j,] = rgamma(n_sample, tests_f_max[i,j+1]*population_f[i,j+1]/100, 1)/pop_active_f[i,j,]
    diag_rate_f_max[i,j,] = rgamma(n_sample, diagnoses_f_max[i,j+1]*population_f[i,j+1]/100000, 1)/pop_active_f[i,j,]
    test_rate_f_min[i,j,] = rgamma(n_sample, tests_f_min[i,j+1]*population_f[i,j+1]/100, 1)/pop_active_f[i,j,]
    diag_rate_f_min[i,j,] = rgamma(n_sample, diagnoses_f_min[i,j+1]*population_f[i,j+1]/100000, 1)/pop_active_f[i,j,]
  }
}


##############################
# function to calculate steady-state A, U and S
##############################

sol_dyn <- list(
  S_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    S <- alpha_AU*alpha_US/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(S)
  },
  
  U_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    U <- alpha_AU*alpha_SU/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(U)
  },
  
  A_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) {
    A <- alpha_SU*alpha_UA/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(A)
  },
  
  prev_fun = function(alpha_UA, alpha_AU, alpha_US, alpha_SU) { #dyn_fun in python code
    prev <- alpha_AU*alpha_US/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA)) +
      alpha_SU*alpha_UA/(alpha_AU*alpha_US + alpha_SU*(alpha_AU + alpha_UA))
    return(prev)
  }
)

##############################
# model the testing and diagnosis rates
##############################

model_test_diag <- list(
  # returns list with two functions
  
  test_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    # function that returns number of tests, given number of asymptomatic cases and uninfecteds
    tsym <- ( ssym + (1 - A - U)*test_sym ) # number of tests
    return(tsym)
  },
  
  diag_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    # function that returns number of diagnoses, given number of asymptomatic cases and uninfecteds
    dsym <- ( A*ssym*true_pos + U*ssym*false_pos + (1 - A - U)*test_sym*true_pos ) # number of diagnoses
    return(dsym)
  }
)

test_diag_fun <- function(parms, p_symp, self_cure, test_sym, true_pos, false_pos, cov=NULL, adpc=NULL){
  # function that returns number of tests and number of diagnoses, given a certain incidence and screening rate (parms),
  # or difference with measured number of tests (cov) and diagnoses (adpc) if cov and adpc not NULL
  inc <- parms[1] # incidence
  scr <- parms[2] # screening rate
  # p_symp = parms[3] # proportion symptomatic
  # self_cure = parms[4] # self-cure rate
  # test_sym = parms[5] # testing rate in symptomatics
  # true_pos = parms[6] # true positive rate
  # false_pos = parms[7] # false positive rate
  
  if (is.null(cov)) cov <- 0
  if (is.null(adpc)) adpc <- 0
  A = sol_dyn$A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
  U = sol_dyn$U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, scr*true_pos + test_sym*true_pos)
  test <- model_test_diag$test_fun(A, U, scr, test_sym, true_pos, false_pos) - cov
  diag <- model_test_diag$diag_fun(A, U, scr, test_sym, true_pos, false_pos) - adpc
  return(c(test=test, diag=diag))
}

##############################
# ...and if symptomatic and asymptomatic diagnoses are observed separately:
##############################

model_test_diag_sa <- list(
  # returns list with two functions
  
  test_sa_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    # function that returns number of tests, given number of asymptomatic cases and uninfecteds
    tsym <- ( (A+U)*ssym + (1 - A - U)*test_sym ) # number of tests
    return(tsym)
  },
  
  diags_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    # function that returns number of symptomatic diagnoses
    dssym <- (1 - A - U)*test_sym*true_pos # number of diagnoses
    return(dsym)
  },
  
  diaga_fun = function(A, U, ssym, test_sym, true_pos, false_pos){
    # function that returns number of asymptomatic diagnoses
    dasym <- ( A*ssym*true_pos + U*ssym*false_pos )
    return(dsym)
  }
)

test_diag_sym_asym_fun <- function(parms, p_symp, self_cure, test_sym, true_pos, false_pos){
  # function that returns number of tests and number of diagnoses, given a certain incidence and screening rate (parms),
  # or difference with measured number of tests (cov) and diagnoses (adpc) if cov and adpc not NULL
  inc <- parms[1] # incidence
  scr <- parms[2] # screening rate
  
  A = sol_dyn$A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
  U = sol_dyn$U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
  test_sa <- model_test_diag_sa$test_sa_fun(A, U, scr, test_sym, true_pos, false_pos)
  diags <- model_test_diag_sa$diags_fun(A, U, scr, test_sym, true_pos, false_pos)
  diaga <- model_test_diag_sa$diaga_fun(A, U, scr, test_sym, true_pos, false_pos)
  return(c(test_sa=test_sa, diags=diags, diaga=diaga))
}

##############################
# use information on symptoms in _prevalent_ (as opposed to incident) infections
##############################

test_diag_prev_symp_fun <- function(parms, p_symp, self_cure, test_sym, true_pos, false_pos){
  inc <- parms[1] # incidence
  scr <- parms[2] # screening rate
  
  A = sol_dyn$A_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
  U = sol_dyn$U_fun(inc*(1-p_symp), self_cure + scr*true_pos, inc*p_symp, self_cure + scr*true_pos + test_sym*true_pos)
  
  test <- model_test_diag$test_fun(A, U, scr, test_sym, true_pos, false_pos)
  diag <- model_test_diag$diag_fun(A, U, scr, test_sym, true_pos, false_pos)
  return(c(test=test, diag=diag, propsym=(1 - A - U)/(1 - U)))
}

######################
# parameters of beta distributions representing the proportion of the population sexually 
# active, by sex and age group
######################

# men, 16-24
alpha_m_16_24 <- get.beta.par(p = c(0.025, 0.975), q = c(0.8023836019, 0.843403825), plot = FALSE)[1]
beta_m_16_24 <- get.beta.par(p = c(0.025, 0.975), q = c(0.8023836019, 0.843403825), plot = FALSE)[2]

# women, 16-24
alpha_f_16_24 <- get.beta.par(p = c(0.025, 0.975), q = c(0.7998634469, 0.837979601), plot = FALSE)[1]
beta_f_16_24 <- get.beta.par(p = c(0.025, 0.975), q = c(0.7998634469, 0.837979601), plot = FALSE)[2]

######################
# sexually-active population:
######################

# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
p_active_m_16_24 = rbeta(n_sample, alpha_m_16_24, beta_m_16_24) # 16-24 yo only
pop_active_m_15_24 = rbinom(n_sample, 3519015, p_active_m_16_24)

p_active_f_16_24 = rbeta(n_sample, alpha_f_16_24, beta_f_16_24) # 16-24 yo only
pop_active_f_15_24 = rbinom(n_sample, 3388842, p_active_f_16_24)

######################
# testing and diagnosis rates, per person per year (Surveillance data on chlamydia testing and diagnosis rates in England in 2012)
# Population, testing and diagnosis data is from http://www.chlamydiascreening.nhs.uk/ps/data.asp (downloaded 17 April 2015).
######################

test_rate_m_15_24 = rgamma(n_sample,566908, 1)/pop_active_m_15_24
diag_rate_m_15_24 = rgamma(n_sample,48387, 1)/pop_active_m_15_24

test_rate_f_15_24 = rgamma(n_sample,1205896, 1)/pop_active_f_15_24
diag_rate_f_15_24 = rgamma(n_sample,88101, 1)/pop_active_f_15_24

######################
# test performance
######################

p_true_pos_m = rbeta(n_sample,32+1, 0+1) # Low Health Technol Assess (2007):
p_false_pos_m = rbeta(n_sample,2+1, 950+1) # Low Health Technol Assess (2007)

p_true_pos_f = rbeta(n_sample,129+1, 12+1) # Low Health Technol Assess (2007): 129 of 141 infected samples tested +ve
p_false_pos_f = rbeta(n_sample,4+1, 2323+1) # Low Health Technol Assess (2007): 4 of 2327 uninfected samples tested +ve

###################################
# Rate of treatment seeking by symptomatic cases
###################################

# 1.2.2 Rate of treatment seeking by symptomatic cases
# (Mercer et al., Sex. Transm. Infect. 83:400-405; 2007).
# Proportion
# Estimate 95% Confidence Interval
# < 1 week 26.7% (14.4, 44.2)%
#   7-13 days 14.4% (6.1, 30.2)%
#   14-27 days 20.8% (13.3, 31.0)%
#   4-6 weeks 16.6% (8.5, 29.9)%
#   > 6 weeks 21.5% (5.5, 56.4)%


# Find beta distributions corresponding to 95% CIs reported in
# Mercer Sex. Transm. Infect. (2007) (see table above).
a = rep(NA, 5)
b = rep(NA, 5)

# < 1 week
a[1] <- get.beta.par(p = c(0.025, 0.975), q = c(0.144, 0.442), plot = FALSE)[1]
b[1] <- get.beta.par(p = c(0.025, 0.975), q = c(0.144, 0.442), plot = FALSE)[2]

# 7-13 days
a[2] <- get.beta.par(p = c(0.025, 0.975), q = c(0.061, 0.302), plot = FALSE)[1]
b[2] <- get.beta.par(p = c(0.025, 0.975), q = c(0.061, 0.302), plot = FALSE)[2]

# 14-27 days
a[3] <- get.beta.par(p = c(0.025, 0.975), q = c(0.133, 0.310), plot = FALSE)[1]
b[3] <- get.beta.par(p = c(0.025, 0.975), q = c(0.133, 0.310), plot = FALSE)[2]

# 28-41 days
a[4] <- get.beta.par(p = c(0.025, 0.975), q = c(0.085, 0.299), plot = FALSE)[1]
b[4] <- get.beta.par(p = c(0.025, 0.975), q = c(0.085, 0.299), plot = FALSE)[2]

# 42 days and over
a[5] <- get.beta.par(p = c(0.025, 0.975), q = c(0.055, 0.564), plot = FALSE)[1]
b[5] <- get.beta.par(p = c(0.025, 0.975), q = c(0.055, 0.564), plot = FALSE)[2]

# Metropolis-Hastings to get a sample for rate of treatment
i = 1

att_symp = rep(NA, n_sample+1000) # testing rate per person per year. Allow 1000 extra samples for burn-in
ll = c(NA,n_sample+1000,5) # log-likelihood
#props = array(NA,dim=c(n_sample+1000,5))
old = 0.04 # starting sample value
new = 0.04 # starting sample value
# simulate probabilities corresponding to data
# proportion expected in each time window
tps = c(0., 7., 14., 28., 42., Inf)
simp_old = exp(-old*tps[1:5]) - exp(-old*tps[2:length(tps)])
simp_new = exp(-new*tps[1:5]) - exp(-new*tps[2:length(tps)])
acc=0.
while(i <= n_sample+1000){# to do samples for p_test_symp
  new = rnorm(1, old, 0.05) # generate a sample from normal distribution
  if (new < 0){
    att_symp[i] = old # reject
    ll[i] = -1e10
  }
  else {
    simp_old = exp(-old*tps[1:5]) - exp(-old*tps[2:length(tps)])
    simp_new = exp(-new*tps[1:5]) - exp(-new*tps[2:length(tps)])
    if (sum(simp_new > 0) != length(tps) - 1){
      att_symp[i] = old # reject
      ll[i] = -1e10
    }
    else{
      # simulate probabilities corresponding to the data
      log_ratio = sum(dbeta((simp_new-0)/1, a, b, log=T) / 1) - sum(dbeta((simp_old-0)/1, a, b, log=T)/1)
      
      if (log(runif(1,0,1)) < log_ratio){
        att_symp[i] = new # accept
        ll[i] = sum(dbeta(simp_new, a, b, log=T))
        old = new
        acc = acc+1
      }
      else{
        att_symp[i] = old # reject
        ll[i] = sum(dbeta(simp_old, a, b))
      }
    }
  }
  #props[i] = simp_old
  i = i+1
}

att_symp = att_symp[1000:length(att_symp)] # remove burn-in samples
ll = ll[1000:length(ll)] # log-likelihood
# acc/(n_sample+1000) # print the proportion of samples accepted
# mean(att_symp)*365.25
# quantile(att_symp, c(.025, .975)) *365.25
att_symp = att_symp*365.25 # convert rate from day^-1 to year^-1

###################################
# Spontaneous clearance
###################################

sc_m  = read.csv('data/chlamydia_two_exponentials_men.csv')[1:n_sample,1]
sc_f  = read.csv('data/chlamydia_two_exponentials_women.csv')[1:n_sample,1]

#########################################################################################################
# Infer the proportion of incident infections which are asymptomatic, by calibrating to the Natsal-3 prevalence estimates.
#########################################################################################################

alpha_prev_m <- get.beta.par(p = c(0.025, 0.975), q = c(0.015, 0.034), plot = FALSE)[1]
beta_prev_m <- get.beta.par(p = c(0.025, 0.975), q = c(0.015, 0.034), plot = FALSE)[2]
prev_m = rbeta(n_sample, alpha_prev_m, beta_prev_m)

alpha_prev_f <- get.beta.par(p = c(0.025, 0.975), q = c(0.022, 0.043), plot = FALSE)[1]
beta_prev_f <- get.beta.par(p = c(0.025, 0.975), q = c(0.022, 0.043), plot = FALSE)[2]
prev_f = rbeta(n_sample, alpha_prev_f, beta_prev_f)

# incidence, screening and proportion of incident infections asymptomatic in men
tmpfun <- function(parms, sc, p_true_pos, p_false_pos, att_symp, mytest_rate=NULL, mydiag_rate=NULL, myprev=NULL){
  # function that returns number of tests& diagnoses and prevalence, given a certain incidence and screening rate,
  # or difference with measured number of tests (mytest_rate) and diagnoses (mydiag_rate) if mytest_rate and mydiag_rate not NULL
  inc <- parms[1]
  scr <- parms[2]
  p_asymp <- parms[3]
  
  if (is.null(mytest_rate)) mytest_rate <- 0
  if (is.null(mydiag_rate)) mydiag_rate <- 0
  if (is.null(myprev)) myprev <- 0
  
  A <- sol_dyn$A_fun(inc*p_asymp, sc + scr*p_true_pos, inc*(1 - p_asymp), scr*p_true_pos + att_symp*p_true_pos)
  U <- sol_dyn$U_fun(inc*p_asymp, sc + scr*p_true_pos, inc*(1 - p_asymp), scr*p_true_pos + att_symp*p_true_pos)
  test <- model_test_diag$test_fun(A, U, scr, att_symp, p_true_pos, p_false_pos) - mytest_rate
  diag <- model_test_diag$diag_fun(A, U, scr, att_symp, p_true_pos, p_false_pos) - mydiag_rate
  
  prev <- sol_dyn$prev_fun(inc*p_asymp,sc + scr*p_true_pos,inc*(1-p_asymp),scr*p_true_pos + att_symp*p_true_pos) - myprev
  
  return(c(test=test, diag=diag, prev=prev))
}

inc_m = rep(0, n_sample)
scr_m = rep(0, n_sample)
p_asymp_m = rep(0, n_sample)

inc_f = rep(0, n_sample)
scr_f = rep(0, n_sample)
p_asymp_f = rep(0, n_sample)

for (i in 1:n_sample){

  # calculate incidence and screening rate such that data on number of tests, diagnoses and prevalence
  # complies with calculated number of tests, diagnoses and prevalence
  solution_m <- multiroot(function(parms) tmpfun(parms, sc_m[i], p_true_pos_m[i], p_false_pos_m[i], att_symp[i], 
                                                 test_rate_m_15_24[i], diag_rate_m_15_24[i], prev_m[i]), 
                          start = c(0.09, 0.25, 0.6))$root
  inc_m[i] <- solution_m[1]
  scr_m[i] <- solution_m[2]
  p_asymp_m[i] <- solution_m[3]
  
}

for (i in 1:n_sample){
    
  solution_f <- multiroot(function(parms) tmpfun(parms, sc_f[i], p_true_pos_f[i], p_false_pos_f[i], att_symp[i], 
                                                 test_rate_f_15_24[i], diag_rate_f_15_24[i], prev_f[i]), 
                          start = c(0.09, 0.25, 0.9))$root
  inc_f[i] <- solution_f[1]
  scr_f[i] <- solution_f[2]
  p_asymp_f[i] <- solution_f[3]
  
}

# The sampled parameter values are now used to infer prevalence in men and women in different age groups, each year.
# max and min refer to the numbers of tests and diagnoses

prev_m_max = array(NA, dim=c(16,2,n_sample))
inc_m_max = array(NA, dim=c(16,2,n_sample))
scr_m_max = array(NA, dim=c(16,2,n_sample))

prev_f_max = array(NA, dim=c(16,2,n_sample))
inc_f_max = array(NA, dim=c(16,2,n_sample))
scr_f_max = array(NA, dim=c(16,2,n_sample))

prev_m_min = array(NA, dim=c(16,2,n_sample))
inc_m_min = array(NA, dim=c(16,2,n_sample))
scr_m_min = array(NA, dim=c(16,2,n_sample))

prev_f_min = array(NA, dim=c(16,2,n_sample))
inc_f_min = array(NA, dim=c(16,2,n_sample))
scr_f_min = array(NA, dim=c(16,2,n_sample))


for (i in 1:16){
  for (j in 1:2){
    for (k in 1:n_sample){
      
      # get incidence and screening rate given data on testing and diagnoses
      ss_m_max <- multiroot(function(parms) test_diag_fun(parms,1-p_asymp_m[k], sc_m[k], att_symp[k], p_true_pos_m[k], p_false_pos_m[k])-c(test_rate_m_max[i,j,k], diag_rate_m_max[i,j,k]),start = c(0.03, 0.44))$root
      inc_m_max[i,j,k] <- ss_m_max[1]
      scr_m_max[i,j,k] <- ss_m_max[2]
  
      # calculate prevalence given the computed incidence and screening rates
      prev_m_max[i,j,k] <- sol_dyn$prev_fun(inc_m_max[i,j,k]*p_asymp_m[k],
                                      sc_m[k] + scr_m_max[i,j,k]*p_true_pos_m[k],
                                      inc_m_max[i,j,k]*(1-p_asymp_m[k]),
                                      scr_m_max[i,j,k]*p_true_pos_m[k] + att_symp[k]*p_true_pos_m[k])
      
      # get incidence and screening rate given data on testing and diagnoses
      ss_f_max <- multiroot(function(parms) test_diag_fun(parms,1-p_asymp_f[k], sc_f[k], att_symp[k], p_true_pos_f[k], p_false_pos_f[k])-c(test_rate_f_max[i,j,k], diag_rate_f_max[i,j,k]),start = c(0.03, 0.44))$root
      inc_f_max[i,j,k] <- ss_f_max[1]
      scr_f_max[i,j,k] <- ss_f_max[2]
  
      # calculate prevalence given the computed incidence and screening rates
      prev_f_max[i,j,k] <- sol_dyn$prev_fun(inc_f_max[i,j,k]*p_asymp_f[k],
                                      sc_f[k] + scr_f_max[i,j,k]*p_true_pos_f[k],
                                      inc_f_max[i,j,k]*(1-p_asymp_f[k]),
                                      scr_f_max[i,j,k]*p_true_pos_f[k] + att_symp[k]*p_true_pos_f[k])
      
      # get incidence and screening rate given data on testing and diagnoses
      ss_m_min <- multiroot(function(parms) test_diag_fun(parms,1-p_asymp_m[k], sc_m[k], att_symp[k], p_true_pos_m[k], p_false_pos_m[k])-c(test_rate_m_min[i,j,k], diag_rate_m_min[i,j,k]),start = c(0.03, 0.44))$root
      inc_m_min[i,j,k] <- ss_m_min[1]
      scr_m_min[i,j,k] <- ss_m_min[2]
  
      # calculate prevalence given the computed incidence and screening rates
      prev_m_min[i,j,k] <- sol_dyn$prev_fun(inc_m_min[i,j,k]*p_asymp_m[k],
                                      sc_m[k] + scr_m_min[i,j,k]*p_true_pos_m[k],
                                      inc_m_min[i,j,k]*(1-p_asymp_m[k]),
                                      scr_m_min[i,j,k]*p_true_pos_m[k] + att_symp[k]*p_true_pos_m[k])
      
      # get incidence and screening rate given data on testing and diagnoses
      ss_f_min <- multiroot(function(parms) test_diag_fun(parms,1-p_asymp_f[k], sc_f[k], att_symp[k], p_true_pos_f[k], p_false_pos_f[k])-c(test_rate_f_min[i,j,k], diag_rate_f_min[i,j,k]),start = c(0.03, 0.44))$root
      inc_f_min[i,j,k] <- ss_f_min[1]
      scr_f_min[i,j,k] <- ss_f_min[2]
  
      # calculate prevalence given the computed incidence and screening rates
      prev_f_min[i,j,k] <- sol_dyn$prev_fun(inc_f_min[i,j,k]*p_asymp_f[k],
                                      sc_f[k] + scr_f_min[i,j,k]*p_true_pos_f[k],
                                      inc_f_min[i,j,k]*(1-p_asymp_f[k]),
                                      scr_f_min[i,j,k]*p_true_pos_f[k] + att_symp[k]*p_true_pos_f[k])
    }
  }
}

prev_m_max_sorted <- aperm (apply (prev_m_max, 1:2, sort) , c(2, 3, 1)) # sort third dimension
prev_f_max_sorted <- aperm (apply (prev_f_max, 1:2, sort) , c(2, 3, 1)) # sort third dimension
prev_m_min_sorted <- aperm (apply (prev_m_min, 1:2, sort) , c(2, 3, 1)) # sort third dimension
prev_f_min_sorted <- aperm (apply (prev_f_min, 1:2, sort) , c(2, 3, 1)) # sort third dimension

# Average duration of infection

durinf_m_max <- prev_m_max/inc_m_max
durinf_m_max_sorted <- aperm (apply (durinf_m_max, 1:2, sort) , c(2, 3, 1)) # sort third dimension
durinf_f_max <- prev_f_max/inc_f_max
durinf_f_max_sorted <- aperm (apply (durinf_f_max, 1:2, sort) , c(2, 3, 1)) # sort third dimension
durinf_m_min <- prev_m_min/inc_m_min
durinf_m_min_sorted <- aperm (apply (durinf_m_min, 1:2, sort) , c(2, 3, 1)) # sort third dimension
durinf_f_min <- prev_f_min/inc_f_min
durinf_f_min_sorted <- aperm (apply (durinf_f_min, 1:2, sort) , c(2, 3, 1)) # sort third dimension

# Per year changes in the transition rates

alpha_UA_f_max <- array(NA, dim=c(16,2,n_sample))
alpha_AU_f_max <- array(NA, dim=c(16,2,n_sample))
alpha_US_f_max <- array(NA, dim=c(16,2,n_sample))
alpha_SU_f_max <- array(NA, dim=c(16,2,n_sample))

alpha_UA_m_max <- array(NA, dim=c(16,2,n_sample))
alpha_AU_m_max <- array(NA, dim=c(16,2,n_sample))
alpha_US_m_max <- array(NA, dim=c(16,2,n_sample))
alpha_SU_m_max <- array(NA, dim=c(16,2,n_sample))

alpha_UA_f_min <- array(NA, dim=c(16,2,n_sample))
alpha_AU_f_min <- array(NA, dim=c(16,2,n_sample))
alpha_US_f_min <- array(NA, dim=c(16,2,n_sample))
alpha_SU_f_min <- array(NA, dim=c(16,2,n_sample))

alpha_UA_m_min <- array(NA, dim=c(16,2,n_sample))
alpha_AU_m_min <- array(NA, dim=c(16,2,n_sample))
alpha_US_m_min <- array(NA, dim=c(16,2,n_sample))
alpha_SU_m_min <- array(NA, dim=c(16,2,n_sample))

# set up a function to simulate system dynamics when perturbed from steady state
for (j in 1:2){
  for (i in 1:16){
    
    # parameters
    alpha_UA_m_max[i,j,] <- inc_m_max[i,j,]*p_asymp_m
    alpha_AU_m_max[i,j,] <- sc_m + scr_m_max[i,j,]*p_true_pos_m
    alpha_US_m_max[i,j,] <- inc_m_max[i,j,]*(1-p_asymp_m)
    alpha_SU_m_max[i,j,] <- scr_m_max[i,j,]*p_true_pos_m + att_symp[1:n_sample]*p_true_pos_m
    
    alpha_UA_m_min[i,j,] <- inc_m_min[i,j,]*p_asymp_m
    alpha_AU_m_min[i,j,] <- sc_m + scr_m_min[i,j,]*p_true_pos_m
    alpha_US_m_min[i,j,] <- inc_m_min[i,j,]*(1-p_asymp_m)
    alpha_SU_m_min[i,j,] <- scr_m_min[i,j,]*p_true_pos_m + att_symp[1:n_sample]*p_true_pos_m
    
    alpha_UA_f_max[i,j,] <- inc_f_max[i,j,]*p_asymp_f
    alpha_AU_f_max[i,j,] <- sc_f + scr_f_max[i,j,]*p_true_pos_f
    alpha_US_f_max[i,j,] <- inc_f_max[i,j,]*(1-p_asymp_f)
    alpha_SU_f_max[i,j,] <- scr_f_max[i,j,]*p_true_pos_f + att_symp[1:n_sample]*p_true_pos_f
    
    alpha_UA_f_min[i,j,] <- inc_f_min[i,j,]*p_asymp_f
    alpha_AU_f_min[i,j,] <- sc_f + scr_f_min[i,j,]*p_true_pos_f
    alpha_US_f_min[i,j,] <- inc_f_min[i,j,]*(1-p_asymp_f)
    alpha_SU_f_min[i,j,] <- scr_f_min[i,j,]*p_true_pos_f + att_symp[1:n_sample]*p_true_pos_f
    
  }
}

alpha_UA_f_max_sorted <- aperm (apply (alpha_UA_f_max, 1:2, sort) , c(2, 3, 1))
alpha_AU_f_max_sorted <- aperm (apply (alpha_AU_f_max, 1:2, sort) , c(2, 3, 1))
alpha_US_f_max_sorted <- aperm (apply (alpha_US_f_max, 1:2, sort) , c(2, 3, 1))
alpha_SU_f_max_sorted <- aperm (apply (alpha_SU_f_max, 1:2, sort) , c(2, 3, 1))
inc_f_max_sorted <- aperm (apply (inc_f_max, 1:2, sort) , c(2, 3, 1))
scr_f_max_sorted <- aperm (apply (scr_f_max, 1:2, sort) , c(2, 3, 1))

alpha_UA_m_max_sorted <- aperm (apply (alpha_UA_m_max, 1:2, sort) , c(2, 3, 1))
alpha_AU_m_max_sorted <- aperm (apply (alpha_AU_m_max, 1:2, sort) , c(2, 3, 1))
alpha_US_m_max_sorted <- aperm (apply (alpha_US_m_max, 1:2, sort) , c(2, 3, 1))
alpha_SU_m_max_sorted <- aperm (apply (alpha_SU_m_max, 1:2, sort) , c(2, 3, 1))
inc_m_max_sorted <- aperm (apply (inc_m_max, 1:2, sort) , c(2, 3, 1))
scr_m_max_sorted <- aperm (apply (scr_m_max, 1:2, sort) , c(2, 3, 1))

alpha_UA_f_min_sorted <- aperm (apply (alpha_UA_f_min, 1:2, sort) , c(2, 3, 1))
alpha_AU_f_min_sorted <- aperm (apply (alpha_AU_f_min, 1:2, sort) , c(2, 3, 1))
alpha_US_f_min_sorted <- aperm (apply (alpha_US_f_min, 1:2, sort) , c(2, 3, 1))
alpha_SU_f_min_sorted <- aperm (apply (alpha_SU_f_min, 1:2, sort) , c(2, 3, 1))
inc_f_min_sorted <- aperm (apply (inc_f_min, 1:2, sort) , c(2, 3, 1))
scr_f_min_sorted <- aperm (apply (scr_f_min, 1:2, sort) , c(2, 3, 1))

alpha_UA_m_min_sorted <- aperm (apply (alpha_UA_m_min, 1:2, sort) , c(2, 3, 1))
alpha_AU_m_min_sorted <- aperm (apply (alpha_AU_m_min, 1:2, sort) , c(2, 3, 1))
alpha_US_m_min_sorted <- aperm (apply (alpha_US_m_min, 1:2, sort) , c(2, 3, 1))
alpha_SU_m_min_sorted <- aperm (apply (alpha_SU_m_min, 1:2, sort) , c(2, 3, 1))
inc_m_min_sorted <- aperm (apply (inc_m_min, 1:2, sort) , c(2, 3, 1))
scr_m_min_sorted <- aperm (apply (scr_m_min, 1:2, sort) , c(2, 3, 1))
