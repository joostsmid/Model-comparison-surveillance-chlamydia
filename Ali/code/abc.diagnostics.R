### R script to give a quick illustration of the posterior credible intervals for key model parameters and ancillary estimates in key age-sex cohorts

### load key modules
source("abc.read.in.data.R") # this will generate a report of 16 "errors" which may be ignored (these errors are due to the missing data in the 2007-2008 test counts)

nyears <- 12

load("output/posterior.dat")
Nsim <- length(epsilon.current)

simulated.population <- simulate.chlamydia(theta)
epsilon <- compute.summary.stats(simulated.population)

#### INFECTIONS

#quartz(width=8,height=6)
pdf(file="~/Desktop/diagnostic.inf.pdf",width=8,height=6)

layout(rbind(c(1,2),c(3,4)))

years <- seq(2001,(2001+nyears-1))

plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.inf.m)),xlab="Year",ylab="# Infections (Total)")

y.matrix <- mock.inf.m[,,2]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"M:15-19",cex=0.5)

y.matrix <- mock.inf.m[,,4]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:25-29",cex=0.5)

y.matrix <- mock.inf.m[,,3]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:20-24",cex=0.5)

y.matrix <- mock.inf.m[,,5]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:30-34",cex=0.5)


plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.inf.f)),xlab="Year",ylab="# Infections (New)")

y.matrix <- mock.inf.f[,,2]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"F:15-19",cex=0.5)

y.matrix <- mock.inf.f[,,4]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:25-29",cex=0.5)

y.matrix <- mock.inf.f[,,3]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:20-24",cex=0.5)

y.matrix <- mock.inf.f[,,5]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:30-34",cex=0.5)



plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.inf.new.m)),xlab="Year",ylab="# Infections (New)")

y.matrix <- mock.inf.new.m[,,2]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"M:15-19",cex=0.5)

y.matrix <- mock.inf.new.m[,,4]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:25-29",cex=0.5)

y.matrix <- mock.inf.new.m[,,3]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:20-24",cex=0.5)

y.matrix <- mock.inf.new.m[,,5]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:30-34",cex=0.5)


plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.inf.new.f)),xlab="Year",ylab="# Infections (Total)")

y.matrix <- mock.inf.new.f[,,2]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"F:15-19",cex=0.5)

y.matrix <- mock.inf.new.f[,,4]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:25-29",cex=0.5)

y.matrix <- mock.inf.new.f[,,3]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:20-24",cex=0.5)

y.matrix <- mock.inf.new.f[,,5]
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:30-34",cex=0.5)

dev.off()


###### PREVALENCE


#quartz(width=8,height=6)
pdf(file="~/Desktop/diagnostic.prev.pdf",width=8,height=6)
layout(rbind(c(1,2),c(3,4)))

years <- seq(2001,(2001+nyears-1))

plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.prev.m*100)),xlab="Year",ylab="Prevalence (%)")

y.matrix <- mock.prev.m[,,2]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"M:15-19",cex=0.5)

y.matrix <- mock.prev.m[,,4]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:25-29",cex=0.5)

y.matrix <- mock.prev.m[,,3]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:20-24",cex=0.5)

y.matrix <- mock.prev.m[,,5]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:30-34",cex=0.5)


plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.prev.f*100)),xlab="Year",ylab="Prevalence (%)")

y.matrix <- mock.prev.f[,,2]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"F:15-19",cex=0.5)

y.matrix <- mock.prev.f[,,4]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:25-29",cex=0.5)

y.matrix <- mock.prev.f[,,3]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:20-24",cex=0.5)

y.matrix <- mock.prev.f[,,5]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:30-34",cex=0.5)




plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.pos.m*100)),xlab="Year",ylab="NNDSS Positivity (%)")

y.matrix <- mock.pos.m[,,2]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"M:15-19",cex=0.5)

y.matrix <- mock.pos.m[,,4]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:25-29",cex=0.5)

y.matrix <- mock.pos.m[,,3]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:20-24",cex=0.5)

y.matrix <- mock.pos.m[,,5]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"M:30-34",cex=0.5)


plot(-1,-1,xlim=c(2001,2001+nyears),ylim=c(0,max(mock.pos.f*100)),xlab="Year",ylab="NNDSS Positivity (%)")

y.matrix <- mock.pos.f[,,2]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=2)
text((2001+nyears-0.25),y.low[nyears],"F:15-19",cex=0.5)

y.matrix <- mock.pos.f[,,4]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=0,border="grey0",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:25-29",cex=0.5)

y.matrix <- mock.pos.f[,,3]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=-45,lwd=1.25,density=10,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:20-24",cex=0.5)

y.matrix <- mock.pos.f[,,5]*100
y.low <- y.high <- numeric(nyears)
for (j in 1:nyears) {
y.low[j] <- quantile(y.matrix[,j],(1-0.95)/2)
y.high[j] <- quantile(y.matrix[,j],1-(1-0.95)/2)}
polygon(c(years,rev(years)),c(y.low,rev(y.high)),angle=45,lwd=1.25,density=15,border="transparent",lty=1)
text((2001+nyears-0.25),y.low[nyears],"F:30-34",cex=0.5)

dev.off()
