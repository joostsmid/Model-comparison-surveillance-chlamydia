### R script to generate plot of positivity surfaces with final model parameters

p.pos <- function(p.inf=0.1,p.screen=0.05,p.asymp=0.92,p.attend.given.symp=0.91,p.test.given.symp=0.95,p.pos=0.97,p.falsepos=0.002) {
    p.pos <- (p.inf*(1-p.asymp)*p.attend.given.symp*p.pos + p.inf*p.asymp*p.screen*p.pos + (1-p.inf)*p.screen*p.falsepos)/(p.inf*(1-p.asymp)*p.attend.given.symp*p.pos + p.inf*p.asymp*p.screen*p.pos + (1-p.inf)*p.screen*p.falsepos + (1-p.inf)*p.screen*(1-p.falsepos))
return(p.pos)
}

p.inf <- seq(0.001,0.15,by=0.001)
p.screen <- seq(0.001,0.25,by=0.001)
p.positive <- matrix(0,nrow=length(p.inf),ncol=length(p.screen))
for (i in 1:length(p.inf)) {for (j in 1:length(p.screen)) {p.positive[i,j] <- p.pos(p.inf=p.inf[i],p.screen=p.screen[j])}}
p.positive.f <- matrix(0,nrow=length(p.inf),ncol=length(p.screen))
for (i in 1:length(p.inf)) {for (j in 1:length(p.screen)) {p.positive.f[i,j] <- p.pos(p.inf=p.inf[i],p.screen=p.screen[j],p.asymp=0.83)}}

setEPS()
postscript("~/Desktop/newfig2.eps",width=6,height=(3.267*2))
#pdf("~/Desktop/newfig2.pdf",width=6,height=(3.267*2))
layout(c(1,2))
par(mai=c(0.45,0.6,0.1,0.1),cex=0.7)

contour(p.screen,p.inf*(1-0.45)*(1-0.09),t(p.positive)*100,levels=c(0.05,0.1,0.15,0.2,0.25,0.3)*100,xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlim=c(0,0.25),ylim=c(0,0.076),col="blue")
#contour(p.screen,p.inf*(1-0.45)*(1-0.09),t(p.positive.f)*100,levels=c(0.05,0.1,0.15,0.2,0.25,0.3)*100,xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlim=c(0,0.25),ylim=c(0,0.076),add=T,col="red")

x <- list()
y <- list()
for (i in c(1,2,3,4,5,8,9,10,11,12,13)) {
p.screenobs <- tested.m[2,i]/(m.15.19+m.20.24)[i]
p.posobs <- colSums(notifications.m[2:3,])[i]/tested.m[2,i]
xx <- yy <- seq(0,0.25,length.out=1000)
for (k in 1:1000) {yy[k] <- p.pos(p.inf=xx[k],p.screen=p.screenobs)}
x[[length(x)+1]] <- p.screenobs
y[[length(y)+1]] <- xx[which.min(abs(yy-p.posobs))]}
x <- as.numeric(x)
y <- as.numeric(y)

points(x,y*(1-0.45)*(1-0.09),pch=19)
text(x,y*(1-0.45)*(1-0.09)+c(0,0,0,0,0,-0.001,0.001,-0.002,0,-0.001,+0.002),c("             2001","             2002","             2003","             2004","             2005","             2008","             2009","             2010","             2011","             2012","             2013"))

axis(1,at=c(0,0.05,0.1,0.15,0.2,0.25),labels=c("0","5","10","15","20","25"),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0,0.05,0.1,0.15,0.2,0.25),labels=c("","","","","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.8,cex=0.9)
axis(2,las=2,at=c(0,0.025,0.05,0.075),labels=c("0","2.5","5","7.5"),tck=-0.0075,hadj=0.5,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.025,0.05,0.075),labels=c("","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=c(0,0.025,0.05,0.075),labels=c("","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Testing Rate (%)",side=1,line=1.65,cex=0.9)

legend("bottomright",c("Males"),lty=1,col=c("blue"),bty='n')
legend("bottom",c("M15-24 including 2001-2011 period of near constant positivity"),pch=19,bty='n')

par(mai=c(0.45,0.6,0.1,0.1),cex=0.7)

#contour(p.screen,p.inf*(1-0.45)*(1-0.09),t(p.positive)*100,levels=c(0.05,0.1,0.15,0.2,0.25,0.3)*100,xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlim=c(0,0.25),ylim=c(0,0.076),col="blue")
contour(p.screen,p.inf*(1-0.45)*(1-0.09),t(p.positive.f)*100,levels=c(0.05,0.1,0.15,0.2,0.25,0.3)*100,xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlim=c(0,0.25),ylim=c(0,0.076),col="red")

#points(c(1.219732,1.477448,1.791748,2.181576,2.334494,3.336556,3.970232,4.542735,5.286047)/100,c(0.036,0.042,0.049,0.057,0.060,0.076,0.085,0.092,0.10)*(1-0.45)*(1-0.09),pch=19)
#text(c(1.219732,1.477448,1.791748,2.181576,2.334494,3.336556,3.970232,4.542735,5.286047)/100,c(0.036,0.042,0.049,0.057,0.060,0.076,0.085,0.092,0.10)*(1-0.45)*(1-0.09),c("             2001","             2002","             2003","             2004","             2005","             2008","             2009","             2010","             2011"))

axis(1,at=c(0,0.05,0.1,0.15,0.2,0.25),labels=c("0","5","10","15","20","25"),tck=-0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
axis(1,at=c(0,0.05,0.1,0.15,0.2,0.25),labels=c("","","","","",""),tck=0.0075,padj=-1.4,cex.axis=1,lwd.ticks=0.5)
mtext("Prevalence (%)",side=2,line=2.8,cex=0.9)
axis(2,las=2,at=c(0,0.025,0.05,0.075),labels=c("0","2.5","5","7.5"),tck=-0.0075,hadj=0.5,cex.axis=1,lwd.ticks=0.5)
axis(2,las=2,at=c(0,0.025,0.05,0.075),labels=c("","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
axis(4,las=2,at=c(0,0.025,0.05,0.075),labels=c("","","",""),tck=0.0075,hadj=0.6,cex.axis=0.8,lwd.ticks=0.5)
mtext("Testing Rate (%)",side=1,line=1.65,cex=0.9)

legend("bottomright",c("Females"),lty=1,col=c("red"),bty='n')
#legend("bottom",c("M15-24 during 2001-2011 period of constant positivity"),pch=19,bty='n')

dev.off()
