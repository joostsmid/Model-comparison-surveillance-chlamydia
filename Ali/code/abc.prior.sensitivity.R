### R script to examine prior-sensitivity in our posteriors estimates for the incidence rate of chlamydia in 15-24 and 25-34 year old men and women

source("abc.read.in.hyperparameters.R")
load("output/posterior2.dat")
load("output/theta.test.50.dat")
library(triangle)

### Varying p_test.|symp.

pi.ratio.lower <- dtriangle(theta[,2],prior.test.given.symp[1],prior.test.given.symp[3],prior.test.given.symp[1])/dtriangle(theta[,2],prior.test.given.symp[1],prior.test.given.symp[3],prior.test.given.symp[2])
pi.ratio.upper <- dtriangle(theta[,2],prior.test.given.symp[1],prior.test.given.symp[3],prior.test.given.symp[3])/dtriangle(theta[,2],prior.test.given.symp[1],prior.test.given.symp[3],prior.test.given.symp[2])

pi.ratio.lower <- pi.ratio.lower/sum(pi.ratio.lower)
pi.ratio.upper <- pi.ratio.upper/sum(pi.ratio.upper)

quantile(sample(mock.incper.m[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.m[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.m[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.m[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.f[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.f[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.f[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.f[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))


quantile(sample(mock.prev.m[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.m[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.m[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.m[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.f[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.f[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.f[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.f[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

### Varying p_pos.|symp.

pi.ratio.lower <- dtriangle(theta[,3],prior.true.pos[1],prior.true.pos[3],prior.true.pos[1])/dtriangle(theta[,3],prior.true.pos[1],prior.true.pos[3],prior.true.pos[2])
pi.ratio.upper <- dtriangle(theta[,3],prior.true.pos[1],prior.true.pos[3],prior.true.pos[3])/dtriangle(theta[,3],prior.true.pos[1],prior.true.pos[3],prior.true.pos[2])

quantile(sample(mock.incper.m[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.m[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.m[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.m[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.f[,12,1]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.f[,12,1]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.incper.f[,12,2]*100,100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.incper.f[,12,2]*100,100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))



quantile(sample(mock.prev.m[,12,1],100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.m[,12,1],100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.m[,12,2],100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.m[,12,2],100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.f[,12,1],100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.f[,12,1],100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))

quantile(sample(mock.prev.f[,12,2],100000,replace=T,prob=pi.ratio.lower),c(0.025,0.5,0.975))
quantile(sample(mock.prev.f[,12,2],100000,replace=T,prob=pi.ratio.upper),c(0.025,0.5,0.975))
