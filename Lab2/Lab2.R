#set working directory 
setwd()

#read in data
Data<- read.csv()
Data_catch <- read.csv()

head(Data)
head(Data_catch)

#Part I
#calculate the proportion at age



#plot age at length


#Calculate catch_at_age_length and catch at age


#plots


#####################
#Part II: Use catch at age data from 3Ps cod to generate a spay plot

#Read in the data


#Plot the age composition data


#generate a 'spay' plot. 
#define functions from Noel

#this will generate the spay plot
bp<-function(x,y,v, scale=3,cohort.id=T,cohort.ps, ...){
  plot(x,y,cex=sqrt(abs(v))*scale, col=ifelse(v<0,'black','grey'), pch=ifelse(v<0,16,16), ...)
  points(x[v>0],y[v>0],cex=sqrt(v[v>0])*scale, col='grey', pch=1, ...)
  
  if(cohort.id){
    lx <- max(x)             
    yval <- unique(y) 
    xval <- unique(x) 
    yval1 <- yval[1:(length(yval)-1)]
    surv.cohort<- lx-yval1
    surv.cohort1 <- substring(as.character(surv.cohort),3,4)
    pos1 <- lx + cohort.ps[1]*diff(range(xval))
    pos2 <- yval1 + cohort.ps[2]*diff(range(yval1))
    text(pos1,pos2,surv.cohort1,cex=0.5,srt=45)
    
    xval1 <- xval[seq(2,length(xval),by=2)]
    ly <- max(y)
    pos1 <- ly + cohort.ps[3]*diff(range(yval))
    pos2 <- xval + cohort.ps[4]*diff(range(xval))
    old.cohort<- xval-ly
    old.cohort1 <- substring(as.character(old.cohort),3,4)
    text(pos2,pos1,old.cohort1,cex=0.5,srt=45)
  }  
}

#this calculates standardized proportions at age per year
spay <- function(tdat){
  stdat <- apply(tdat,2,sum,na.rm=T)
  ny <- length(stdat)
  na <- nrow(tdat)
  mstdat <- matrix(stdat,nrow=na,ncol=ny,byrow=T)
  ptdat <- tdat/mstdat
  ptdat.dev <- ptdat - matrix(apply(ptdat,1,mean,na.rm=T),nrow=na,ncol=ny,byrow=F)
  ptdat.std <- sqrt(apply(ptdat.dev^2,1,mean,na.rm=T))
  ztdat <- ptdat.dev/matrix(ptdat.std,nrow=na,ncol=ny,byrow=F)
  return(ztdat)}

scale=1.3

#apply to our data 
sdat <- t(spay(t(Cm))) 
syear <- rep(data[,1],n.age)
sage <- rep(3:14,each=n.year)
vspay <- as.vector(unlist(sdat))

cohort.ps=c(0.02,0.01,0.03,0.01)
y.adj = 0.2
ylim=range(sage);
ylim[2]=ylim[2]+y.adj
xlim=range(syear);
xlim[2]=xlim[2]+0.2

#produce spay plot
par(mfcol=c(1,1),mar=c(3,3,2,2),las=1,mgp=c(2,1,0),oma=c(1,1,1,1))
bp(syear,sage,vspay,cohort.ps=cohort.ps, ylim=ylim, xlim=xlim, 
   las=1, xlab='', ylab='', main="", scale=scale)
points(syear,sage,cex=sqrt(abs(vspay))*scale, col='black', pch=1)
mtext(side=3,outer=F,line=.1,c('3Ps Cod C@A SPAY'),cex=1) 
mtext(side=2,outer=F,line=2,'Age',cex=1)  
mtext(side=1,outer=F,line=2,'Year',cex=1)


#replicate the spay plot using ggplot 
library(reshape2)
library(tidyverse)
