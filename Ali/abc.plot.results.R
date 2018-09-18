### R script for generating the key chlamydia plots of the H. Ali et al. paper: Figs 1, 2[a-d] & 3[a-b]

### this script requires manual editing by someone who understands R in order to plot results beyond 2012 [ignoring our current extrapolation to 2013]

# Output folder - assuming correct working directory has been set
outputFolder <- file.path(getwd(),"output","figures")

# Specify type of image file for saved plots - pdf or eps
plotpdf <- TRUE

# Specify doing all plots 
allPlots <- TRUE

myquantile <- function(x) quantile(x, probs = c(0.025, .5, 0.975), na.rm = T)

source("code/load.library.R")

if (allPlots) {

  # load(file.path("C:/Users/smid/Desktop/chlamydia UK/output",paste0("posterior_",dataset,".dat")))
  load(file.path(getwd(),"output",paste0("posterior_",dataset,".dat")))

  source("code/abc.read.in.data.R")
  
  setEPS()
  if (plotpdf) {
    pdf(file.path(outputFolder,paste0("fig3_",dataset,".pdf")),width=7,height=10)
  } else {
    postscript(file.path(outputFolder,paste0("fig3_",dataset,".eps")),width=6,height=3.267*3)
  }
  
  layout(cbind(c(1,2,3)))
  
  ## Fig 3[a]: Notifications by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,5*10^4),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.f[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.f[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=1,col="grey79",border="transparent",density=-1)

  y.matrix <- mock.not.m[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.not.m[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of notifications",side=2,line=3.4,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("           0",""," 10,000",""," 20,000",""," 30,000",""," 40,000",""," 50,000"),tck=-0.0075,hadj=0.77,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2000:(2000+nyears-1),notifications.f[1,],pch=19,cex=0.8)
  points(2000:(2000+nyears-1),notifications.f[2,],pch=21,cex=0.9)
  points(2000:(2000+nyears-1),notifications.m[1,],pch=15,cex=0.8)
  points(2000:(2000+nyears-1),notifications.m[2,],pch=22,cex=0.9)
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  ## Fig 3[b]: Test count by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,8*10^5),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.f[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.f[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=1,col="grey79",border="transparent",density=-1)
  
  y.matrix <- mock.test.m[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.test.m[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of tests",side=2,line=3.7,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("             0",""," 100,000",""," 200,000",""," 300,000",""," 400,000","","500,000","","600,000"),tck=-0.0075,hadj=0.8,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2000:(2000+nyears-1),tested.f[1,],pch=19,cex=0.8)
  points(2000:(2000+nyears-1),tested.f[2,],pch=21,cex=0.9)
  points(2000:(2000+nyears-1),tested.m[1,],pch=15,cex=0.8)
  points(2000:(2000+nyears-1),tested.m[2,],pch=22,cex=0.9)
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  ## Fig 3[c]: Prevalence by sex & age-group: 15-24 & 25-29
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  xticks <- seq(2000,(2000+nyears-1),by=2)
  finalYear <- (2000+nyears-1)
  axis(1,at=xticks,tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=xticks,labels=rep("",length(xticks)),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
  axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  box()
  
  dev.off()
  
  ##########################################################################################################
  ### Fig. 3a, as in Lewis-White
  ##########################################################################################################
  
  setEPS()
  if (plotpdf) {
    pdf(file.path(outputFolder,paste0("fig3a_",dataset,".pdf")),width=7,height=10)
  } else {
    postscript(file.path(outputFolder,paste0("fig3_",dataset,".eps")),width=6,height=3.267*3)
  }
  
  layout(cbind(c(1,2,3)))
  
  ## Fig 3[a]: Notifications by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,35),xlab="",ylab="",xaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.not.f[,,1]
  y.matrix <- t(t(mock.not.f[,,1]) / f.15.19) * 1000
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.not.f[,,2]
  y.matrix <- t(t(mock.not.f[,,2]) / f.20.24) * 1000
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=1,col="grey79",border="transparent",density=-1)
  
  #y.matrix <- mock.not.m[,,1]
  y.matrix <- t(t(mock.not.m[,,1]) / m.15.19) * 1000
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.not.m[,,2]
  y.matrix <- t(t(mock.not.m[,,2]) / m.20.24) * 1000
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of notifications",side=2,line=3.4,cex=0.9)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^4,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2000:(2000+nyears-1),notifications.f[1,]/f.15.19*1000,pch=19,cex=0.8)
  points(2000:(2000+nyears-1),notifications.f[2,]/f.20.24*1000,pch=21,cex=0.9)
  points(2000:(2000+nyears-1),notifications.m[1,]/m.15.19*1000,pch=15,cex=0.8)
  points(2000:(2000+nyears-1),notifications.m[2,]/m.20.24*1000,pch=22,cex=0.9)
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  ## Fig 3[b]: Test count by sex & age-group: 15-24 & 25+
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,45),xlab="",ylab="",xaxt='n')#,yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.test.f[,,1]
  y.matrix <- t(t(mock.test.f[,,1]) / f.15.19) * 100
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.test.f[,,2]
  y.matrix <- t(t(mock.test.f[,,2]) / f.20.24) * 100
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=1,col="grey79",border="transparent",density=-1)
  
  #y.matrix <- mock.test.m[,,1]
  y.matrix <- t(t(mock.test.m[,,1]) / m.15.19) * 100
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.test.m[,,2]
  y.matrix <- t(t(mock.test.m[,,2]) / m.20.24) * 100
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Number of tests",side=2,line=3.7,cex=0.9)
  # axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("             0",""," 100,000",""," 200,000",""," 300,000",""," 400,000","","500,000","","600,000"),tck=-0.0075,hadj=0.8,cex.axis=1,lwd.ticks=0.5)
  # axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  # axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)*10^5,labels=c("            ","","","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  points(2000:(2000+nyears-1),tested.f[1,]/f.15.19*100,pch=19,cex=0.8)
  points(2000:(2000+nyears-1),tested.f[2,]/f.20.24*100,pch=21,cex=0.9)
  points(2000:(2000+nyears-1),tested.m[1,]/m.15.19*100,pch=15,cex=0.8)
  points(2000:(2000+nyears-1),tested.m[2,]/m.20.24*100,pch=22,cex=0.9)
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  ## Fig 3[c]: Prevalence by sex & age-group: 15-24 & 25-29
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,0.075),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.f[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.prev.m[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  xticks <- seq(2000,(2000+nyears-1),by=2)
  finalYear <- (2000+nyears-1)
  axis(1,at=xticks,tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=xticks,labels=rep("",length(xticks)),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Prevalence (%)",side=2,line=2.5,cex=0.9)
  axis(2,las=2,at=seq(0,0.075,by=0.0125),labels=c("     0","1.25","  2.5","3.75","     5","6.25","  7.5"),tck=-0.0075,hadj=0.55,cex.axis=1,lwd.ticks=0.5)
  axis(4,las=2,at=seq(0,0.075,by=0.0125),labels=c("  ","  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  box()
  
  dev.off()
  
  ##########################################################################################################
  ### Fig. 4
  ##########################################################################################################
  
  # Create plot
  setEPS()
  if (plotpdf) {
    pdf(file.path(outputFolder,paste0("fig4_",dataset,".pdf")),width=7,height=10)
  } else {
    postscript(file.path(outputFolder,paste0("fig4_",dataset,".eps")),width=6,height=3.267*2)
  }
  
  layout((c(1,2)))
  par(mai=c(0.45,0.90,0.10,0.10),cex=0.7)
  
  maxy <- 2
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)),ylim=c(0,maxy*10^5),xlab="",ylab="",xaxt='n',yaxt='n')
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.f[,,1] 
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.f[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.m[,,1]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.inc.m[,,2]
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Estimated annual number of incident",side=2,line=4.8,cex=0.9)
  mtext("chlamydia cases",side=2,line=3.55,cex=0.9)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5)*10^5,labels=c("             0"," 50,000"," 100,000","150,000"," 200,000",""),tck=-0.0075,hadj=0.8,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^5,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)*10^5,labels=c("            ","","","","","","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  ## Fig 4[b]: Incidence percentage by sex & age-group
  par(mai=c(0.45,0.7,0.10,0.15),cex=0.7)
  
  plot(-100,-100,xlim=c(2000,(2000+nyears-1)+1),ylim=c(0,0.12),xlab="",ylab="",xaxt='n',yaxt='n')
  
  med.incper.f1 <- rep(NA, nyears)
  med.incper.f2 <- rep(NA, nyears)
  med.incper.m1 <- rep(NA, nyears)
  med.incper.m2 <- rep(NA, nyears)
  
  
  for (i in 1:nyears) {
    
    med.incper.f1[i] <- myquantile(mock.inc.f[,i,1])[2]/f.15.19[i]
    med.incper.f2[i] <- myquantile(mock.inc.f[,i,2])[2]/f.20.24[i]
    med.incper.m1[i] <- myquantile(mock.inc.m[,i,1])[2]/m.15.19[i]
    med.incper.m2[i] <- myquantile(mock.inc.m[,i,2])[2]/m.20.24[i]
  }
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.incper.f[,,1]
  y.matrix <- t(t(mock.inc.f[,,1]) / f.15.19)
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey89",border="transparent",density=-1)
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.incper.f[,,2]
  y.matrix <- t(t(mock.inc.f[,,2]) / f.20.24)
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=90,lwd=0.5,col="grey79",border="transparent",density=-1)
  
  #y.low <- y.high <- numeric(nyears)
  y.matrix <- mock.incper.m[,,1]
  y.matrix <- t(t(mock.inc.m[,,1]) / m.15.19)
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=1.25,col="grey15",density=10,border="transparent")
  
  y.low <- y.high <- numeric(nyears)
  #y.matrix <- mock.incper.m[,,2]
  y.matrix <- t(t(mock.inc.m[,,2]) / m.20.24)
  for (j in 1:nyears) {
    y.low[j] <- myquantile(y.matrix[,j])[1]
    y.high[j] <- myquantile(y.matrix[,j])[3]
  }
  polygon(c(seq(2000,(2000+nyears-1),by=1),rev(seq(2000,(2000+nyears-1),by=1))),c(y.low,rev(y.high)),
          angle=-45,lwd=0.5,col="grey15",density=30,border="transparent")
  
  points(2000:(2000+nyears-1),med.incper.f1,pch=19,cex=0.8)
  points(2000:(2000+nyears-1),med.incper.f2,pch=21,cex=0.9)
  points(2000:(2000+nyears-1),med.incper.m1,pch=15,cex=0.8)
  points(2000:(2000+nyears-1),med.incper.m2,pch=22,cex=0.9)
  
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
  axis(1,at=c(2000,2002,2004,2006,2008,2010,2012,2014),labels=c("","","","","","","",""),tck=0.0075,padj=-1.7,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Annual incidence (%)",side=2,line=2.4,cex=0.9)
  axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("     0","  2.5","     5","  7.5"," 10.0"," 12.5"),tck=-0.0075,hadj=0.6,cex.axis=1,lwd.ticks=0.5)
  axis(2,las=2,at=seq(0,0.125,by=0.025),labels=c("    ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  axis(4,las=2,at=seq(0,0.125,by=0.025),labels=c("  ","  ","   ","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
  mtext("Year",side=1,line=1.65,cex=0.9)
  
  
  legend("bottomright",c("F 15-19","F 20-24","M 15-19","F 20-24"),pch=c(19,21,15,22),bty='n',cex=1,pt.cex=c(0.8,0.8,0.9,0.9))
  
  dev.off()
  
}

##########################################################################################################
### Fig 5: Time Independent Parameter Constraints
##########################################################################################################

source("code/abc.read.in.hyperparameters.R")

setEPS()
if (plotpdf) {
  pdf(file.path(outputFolder,paste0("fig5_",dataset,".pdf")),width=6,height=6)
} else {
  postscript(file.path(outputFolder,paste0("fig5_",dataset,".eps")),width=6,height=6)
}


par(mai=c(0.45,0.6,0.10,0.15),cex=0.7)

layout(rbind(c(1,2,3),c(4,5,6),c(7,8,9)))

### Asymp.m rbeta(Nsim,400,40) # asymp.m [0.881,0.910,0.934]
posterior <- density(theta[,6],from=0,to=1,n=512,bw="SJ")

plot(-100,-100,xlim=c(0.86,0.95),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.86,0.95,length.out=100)
lines(xx,dbeta(xx,400,40)/max(dbeta(xx,400,40)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.88,0.91,0.94),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.88,0.91,0.94),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[asymp*plain(". [M]")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### Asymp.f rbeta(Nsim,400,75) # asymp.f [0.808,0.843,0.873]
posterior <- density(theta[,7],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.79,0.89),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.79,0.89,length.out=100)
lines(xx,dbeta(xx,400,75)/max(dbeta(xx,400,75)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.80,0.83,0.86),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.80,0.83,0.86),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[asymp*plain(". [F]")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### Attend.symp # attend|symp [0.863,0.905,0.939] # 4000 420 [0.896,0.905,0.913]
posterior <- density(theta[,1],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.85,0.95),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.85,0.95,length.out=100)
lines(xx,dbeta(xx,4000,420)/max(dbeta(xx,4000,420)),lwd=1,col="grey79")
lines(xx,dbeta(xx,207,22)/max(dbeta(xx,207,22)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,1],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,4000,320)/max(dbeta(xx,4000,320)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,3200,420)/max(dbeta(xx,3200,420)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.86,0.9,0.94),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.86,0.9,0.94),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[attend*plain("|symp.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.test.given.symp # rbeta(Nsim,83,4) # test|symp [0.901,0.957,0.987]

posterior <- density(theta[,2],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.88,1),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.88,1,length.out=100)
lines(xx,dbeta(xx,83,4)/max(dbeta(xx,83,4)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

lines(xx,dbeta(xx,500,10)/max(dbeta(xx,500,10)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,1000,75)/max(dbeta(xx,1000,75)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.91,0.95,0.99),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.91,0.95,0.99),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[test*plain("|symp.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

## true.pos # rbeta(Nsim,110,7) # true.pos [0.891,0.943,0.975]

posterior <- density(theta[,3],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.87,0.99),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.87,0.99,length.out=100)
lines(xx,dbeta(xx,110,7)/max(dbeta(xx,83,4)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

lines(xx,dbeta(xx,410,10)/max(dbeta(xx,410,10)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,810,70)/max(dbeta(xx,810,70)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.88,0.93,0.98),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.88,0.93,0.98),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[pos.]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### false.pos rbeta(Nsim,2,250) # false.pos [0.001,0.007,0.022]

posterior <- density(theta[,4],from=0,to=1,n=2048,kernel="gaussian")

plot(-100,-100,xlim=c(0,0.025),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0,0.025,length.out=100)
lines(xx,dbeta(xx,2,250)/max(dbeta(xx,2,250)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0,0.01,0.02),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0,0.01,0.02),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[false*plain(" pos.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### rep rbeta(Nsim,150,1) # p.rep [0.976,0.995,0.999] # 1500 8 [0.990,0.995,0.998]
posterior <- density(theta[,5],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.98,1),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.98,1,length.out=100)
lines(xx,dbeta(xx,1500,8)/max(dbeta(xx,1500,8)),lwd=1,col="grey79")
lines(xx,dbeta(xx,150,1)/max(dbeta(xx,150,1)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,5],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,800,1)/max(dbeta(xx,800,4)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,3000,30)/max(dbeta(xx,3000,30)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.98,0.99,1),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.98,0.99,1),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[rep.]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.cured.after.year # rbeta(Nsim,22,26) # p.cured.after.year [0.321,0.458,0.599]

posterior <- density(theta[,8],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.25,0.65),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.25,0.65,length.out=100)
lines(xx,dbeta(xx,22,26)/max(dbeta(xx,22,26)),lwd=1,col="grey79")

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")

axis(1,at=c(0.30,0.45,0.6),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.3,0.45,0.6),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[self*plain(" cure, n.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

### prior.cured.after.year # rbeta(Nsim,29,290) # p.background.antibiotic [0.062,0.090,0.125] # 1050 10000 [0.090,0.095,0.101]
posterior <- density(theta[,9],from=0,to=1,n=2048)

plot(-100,-100,xlim=c(0.05,0.14),ylim=c(0,1.2),xlab="",ylab="",xaxt='n',yaxt='n',xaxs='i',yaxs='i')

xx <- seq(0.05,0.14,length.out=100)
lines(xx,dbeta(xx,1050,10000)/max(dbeta(xx,1050,10000)),lwd=1,col="grey79")
lines(xx,dbeta(xx,29,290)/max(dbeta(xx,29,290)),lwd=1,col="grey79",lty=2)

lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15")
posterior <- density(theta[,9],from=0,to=1)
lines(posterior$x,posterior$y/max(posterior$y),lwd=1.25,col="grey15",lty=2)

lines(xx,dbeta(xx,270,2000)/max(dbeta(xx,270,2000)),lwd=1,col="red",lty=3)
lines(xx,dbeta(xx,150,2000)/max(dbeta(xx,150,2000)),lwd=1,col="blue",lty=3)

axis(1,at=c(0.05,0.09,0.13),tck=-0.0075,padj=-1.3,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0.05,0.09,0.13),labels=c("","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Relative density",side=2,line=2.4,cex=0.9)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),tck=-0.0075,hadj=0.47,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("","","","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext(expression(p[self*plain(" cure, a.")]),side=1,line=1.65,cex=0.9)

legend("topright",c("Prior"),lwd=c(1,1.25),col=c("grey79","grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)
legend("topleft",c("Post."),lwd=c(1.25),col=c("grey15"),ncol=1,bty='n',cex=1,x.intersp=0.8)

dev.off()

