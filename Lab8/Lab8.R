# Lab 8
# hake
library("nlstools")
library(tidyverse)

#for use in illustrative examples
# q = 0.0004360
# r = 0.379
# K = 4000
# catch = hake.data$catch

fname = c('hake_data.txt')
hake.data = read.table(fname,header=T)
hake.data$log.index = log(hake.data$index)

spm = function(r,K,catch){
  n=length(catch)
  B=rep(NA,n);
  B[1]=K
  for (i in 2:n){
    B[i]=B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catch[i-1]
  }
  B[B<0]=1e-10
  return(B)
}

spm_fit = function(logq,logr,logK,catch){
  B = spm(exp(logr),exp(logK),catch)
  return(data.frame(B=B,EI=exp(logq)*B))
}

spm.fit <- nls(log.index ~ log(spm_fit(logq,logr,logK,catch)[,2]), data = hake.data,
               start = list(logq = -5, logr=log(0.3), logK = log(4000)))
summary(spm.fit)
spm.profile = confint(spm.fit)
print(spm.profile)

parm.est = exp(coef(spm.fit))
names(parm.est)=c("q","r","K")
rownames(spm.profile) =  c("q","r","K")
print(parm.est)
round(parm.est, 4)
print(exp(spm.profile))

nr <- nlsResiduals(spm.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

n = length(hake.data$index)
spm.predict = exp(predict(spm.fit))

par(mfrow=c(2,1),mar=c(0,3,0.5,1),oma=c(3,0,0,0),mgp=c(2,1,0))
plot(hake.data$Year,hake.data$index,xlab='Year',
     ylab='Index',type='p',lwd=2,xaxt='n')
lines(hake.data$Year,spm.predict,lwd=2,col='red')                         

plot(hake.data$Year,nr$resi2[,2],xlab='Year',
     ylab='Std residuals',type='p',lwd=2)
abline(h=0,lty=2)

parm.est = coef(spm.fit)
pop.est = spm_fit(parm.est[1],parm.est[2],parm.est[3],hake.data$catch)

par(mar=c(3.5,3.5,1,1),mgp=c(2.5,1,0))
plot(hake.data$Year,pop.est[,1],xlab='Year',
     ylab='Biomass',type='b',lwd=2)
Bmsy = exp(parm.est[3])/2
abline(h=Bmsy,lty=2,lwd=2)
text(hake.data$Year[1]+1,Bmsy+0.03*diff(range(pop.est[,1])),"Bmsy")

#####################################;

spm = function(r,K,catch){
  n=length(catch)
  B=rep(NA,n);
  B[1]=2*K
  for (i in (2:n)){
    B[i]=B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catch[i-1]
  }
  B[B<0]=1e-10
  return(B)
}

spm.fit1 <- nls(log.index ~ log(spm_fit(logq,logr,logK,catch)[,2]), data = hake.data,
                start = list(logq = -5, logr=log(0.3), logK = log(4000)))
summary(spm.fit1)
spm.profile = confint(spm.fit1)

print(spm.profile)

parm.est = exp(coef(spm.fit1))
names(parm.est)=c("q","r","K")

rownames(spm.profile) =  c("q","r","K")
print(parm.est)
print(exp(spm.profile))

parm.est = coef(spm.fit)
pop.est = spm_fit(parm.est[1],parm.est[2],parm.est[3],hake.data$catch)

par(mar=c(3.5,3.5,1,1),mgp=c(2.5,1,0))
plot(hake.data$Year,pop.est[,1],xlab='Year',
     ylab='Biomass',type='b',lwd=2)

Bmsy = exp(parm.est[3])/2
abline(h=Bmsy,lty=2,lwd=2)
text(hake.data$Year[1]+1,Bmsy+0.03*diff(range(pop.est[,1])),"Bmsy")


nr <- nlsResiduals(spm.fit1)
par(mar=c(4,4,3,2),mgp=c(3,1,0))
plot(nr, which = 0)

n = length(hake.data$index)
spm.predict = exp(predict(spm.fit1))

par(mfrow=c(2,1),mar=c(0,3,0.5,1),oma=c(3,0,0,0),mgp=c(2,1,0))
plot(hake.data$Year,hake.data$index,xlab='Year',
     ylab='Index',type='p',lwd=2,xaxt='n')
lines(hake.data$Year,spm.predict,lwd=2,col='red')

plot(hake.data$Year,nr$resi2[,2],xlab='Year',
     ylab='Std residuals',type='p',lwd=2)
abline(h=0,lty=2)

test.nlsResiduals(nr)

## stock status relative to Bmsy
boo <- nlsBoot(spm.fit1) 
par(oma=c(0,0,0,0),mar=c(2,2,1,1),mgp=c(2.5,1,0))
plot(boo)

par(mfrow=c(2,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
hist(boo$coef[,2],xlab="r",main='')
hist(boo$coef[,3],xlab="K",main='')

options(digits=3)
summary(boo)

spm.wrap = function(parm,catch){
  r=exp(parm[1])
  K=exp(parm[2])
  B = spm(r,K,catch)
  stock.status =5*B[length(B)]/K
  return(stock.status)
}

boo.status = apply(boo$coef[,2:3],1,spm.wrap,
                   catch=hake.data$catch);

par(mar=c(3.5,3.5,1,1),mgp=c(2.5,1,0))
hist(boo.status,xlab="B1988/40%Bmsy",main='')

options(digits=3)
quantile(boo.status,probs=c(0.5,0.025,0.975))

#-------------------------------------------
# redfish
library(R2HTML)
library("nlstools")

Cdat <- read.table("catch.txt",header=TRUE, fill = FALSE)
Cdat$imap = 1:nrow(Cdat); ## year indicator
Cdat$catch = Cdat$catch/1000 ## units in Kilo tonnes
head(Cdat)

Idat <- read.table("indices.txt",header=TRUE, fill = FALSE)
Idat$iq = as.numeric(as.factor(Idat$name))
Idat$iyear = Cdat$imap[match(Idat$year,Cdat$year)] 
Idat$log_index = log(Idat$index)
head(Idat)

spm = function(r,K,catch){
  n=length(catch)
  B=rep(NA,n);
  B[1]=K/2
  for (i in 2:n){
    B[i]=B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catch[i-1]
  }
  B[B<0]=1e-10
  return(B)
}

spm_fit = function(logq,logr,logK,log_index,iq,iyear){
  B = spm(exp(logr),exp(logK),Cdat$catch)
  log_EIndex = logq[iq] + log(B[iyear])
  resid = log_index - log_EIndex
  return(cbind(log_EIndex,resid))
}

n.index = length(unique(Idat$iq))

start.logq = aggregate(Idat$log_index,list(name=Idat$iq),mean)$x - log(400)
sparm = list(logq=start.logq,logr=log(0.1),logK=log(800))

spm.fit <- nls(log_index ~ spm_fit(logq,logr,logK,log_index,iq,iyear)[,1],
               data = Idat,
               start = sparm,
               algorithm="port",
               control=list(maxiter=1000),
               lower=c(rep(-Inf,n.index),log(0.01),log(100)),
               upper=c(rep(Inf,n.index),log(0.3),log(2000)))

summary(spm.fit)

spm.profile = confint.default(spm.fit)
print(spm.profile)

parm.est = exp(coef(spm.fit))
index.name = aggregate(Idat$name,list(iq=Idat$iq),unique)$x
parm.name = c(paste("q",index.name),"r","K")
out = cbind(parm.est,exp(spm.profile))
rownames(out) =  parm.name

print(out)
round(out,3)

########residuals
nr <- nlsResiduals(spm.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

parm.est = coef(spm.fit)
pop.est = spm(exp(parm.est[9]),exp(parm.est[10]),Cdat$catch) 
Bmsy = exp(parm.est[10])/2

par(mar=c(3.5,3.5,1,1),mgp=c(2.5,1,0))
plot(Cdat$year,pop.est/Bmsy,xlab='Year',
     ylab='B/Bmsy',type='b',lwd=2)
abline(h=1,lty=2,lwd=2)

################# CPUE plot ############################;
require("lattice")    
require("xtable")
index.dat <- Idat
rep = data.frame(Eindex=exp(nr$resi3[,1]),std_resid=nr$resi2[,2])


uindex = unique(index.dat$name)
posv = c('topright',rep("topleft",7))

for(i in 1: length(uindex)){
  
  u = uindex[i]
  pos = posv[i]
  
  scale=1000
  index.name = "Index (000s)"
  if(u=='CPUE'){scale=1;index.name="Index"}
  
  ind = index.dat$name==u
  year = index.dat$year[ind]
  index = index.dat$index[ind]/scale
  Eindex = rep$Eindex[ind]/scale
  std.resid = rep$std_resid[ind]
  
  jpeg(file=paste('fit_',u,'.jpeg',sep=''),width=3,height=4,units='in',res=300)
  
  par(mfcol=c(3,1),oma=c(1,1,0.5,1),mar=c(3,3,1,0),las=1,cex.lab=0.7)
  
  ylim=range(index,Eindex)
  xlim=range(index.dat$year)
  
  xi = Cdat$year
  yi=rep(NA,length(xi));pyi=yi;pres=yi
  ind1 = xi%in%year
  yi[ind1]=index 
  pyi[ind1]=Eindex 
  pres[ind1]=std.resid
  
  plot(xi,yi,type='b',lwd=2,pch=19,xlab='',ylab='',ylim=ylim,xlim=xlim)
  lines(xi,pyi,type='l',lwd=2,col='red')   
  points(xi,pyi,col='red')
  mtext(side=1,line=2,outer=F,"Year")            
  mtext(side=2,line=2.5,outer=F,index.name,las=0)
  legend(pos,bty='n',lwd=2,pch=c(19,NA),lty=1,c("Observed","Predicted"),col=c('black','red'))
  
  plot(xi,pres,type='b',xlab='',ylab='',xlim=xlim)
  abline(h=0,lty=2)
  mtext(side=1,line=2,outer=F,"Year")            
  mtext(side=2,line=2,outer=F,"Std resid",las=0)
  
  plot(pyi,pres,type='p',xlab='',ylab='')
  mtext(side=1,line=2,outer=F,"Predicted")            
  mtext(side=2,line=2,outer=F,"Std resid",las=0)
  
  mtext(side=3,line=-1,outer=T,u,las=0)
  
  dev.off() 
}

tble=sqrt(tapply(nr$resi1[,2],index.dat$name,var))

par(mar=c(4,6,1,1),mgp=c(2,1,0))
barplot(tble,horiz = TRUE,las=1,xlab='Residual Std')

##############  biomass and harvest plot ################;
biomass = pop.est;
exploit = Cdat$catch/biomass;
tmb.data=Cdat
Hmsy = exp(parm.est[9])/2

jpeg(file="pop.jpeg",width=4,height=5,units='in',res=1200)

par(mfcol=c(2,1),oma=c(1,1,0,0.5),mar=c(2,3,1,0),cex.axis=1,las=1)
xlim <- c(min(tmb.data$year),max(tmb.data$year))
ylim <- range(biomass)
n=length(tmb.data$year)

ytext <- c("Exploitable biomass (Kt)")
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)
lines(tmb.data$year,biomass,type='l',lty=1,lwd=2)  
abline(h=Bmsy,lwd=2,col='red')    
abline(h=0.3*Bmsy,lwd=1,col='red')    
mtext(side=2,line=3,outer=F,ytext,cex=1,las=0)
p1 = Bmsy + 0.05*diff(range(ylim))
text(tmb.data$year[1]+2,p1,"Bmsy",cex=0.75)   
p1 = 0.3*Bmsy + 0.05*diff(range(ylim))
text(tmb.data$year[1]+2,p1,"Blim",cex=0.75)


ylim <- c(min(exploit),min(400,max(exploit)))
ylim[2]=0.5

ytext <- c("Exploitation rate") 
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)
lines(tmb.data$year,exploit,type='l',lty=1,lwd=2)  
abline(h=Hmsy,lwd=2,col='red')
mtext(side=2,line=3,outer=F,ytext,cex=1,las=0)  
p1 = Hmsy + 0.05*diff(range(ylim))
text(tmb.data$year[n]-6,p1,"Flim=Hmsy",cex=0.75)
mtext(side=1,line=2,outer=F,"Year",cex=1.2)

dev.off()


##############  biomass and harvest wrt MSY plot ################;
rbm=biomass/Bmsy
rF = exploit/Hmsy

jpeg(file="rpop.jpeg",width=4,height=5,units='in',res=1200)

par(mfcol=c(2,1),oma=c(1,0.5,0,0.5),mar=c(2,3,1,0),cex.axis=1,las=1)
xlim <- c(min(tmb.data$year),max(tmb.data$year))
ylim <- c(min(rbm),max(rbm))

n=length(tmb.data$year)

ytext <- c("Biomass/Bmsy")
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)
lines(tmb.data$year,rbm,type='l',lty=1,lwd=2)     
abline(h=1,lwd=2,col='red')   
abline(h=0.3,lwd=1,col='red')       
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)

ylim <- c(min(rF),min(400,max(rF)))
ytext <- c("H/Hmsy") 
plot(xlim,ylim,xlab='',ylab='',type='n',cex.lab=1,lwd=2)
lines(tmb.data$year,rF,type='l',lty=1,lwd=2)  
abline(h=1,lwd=2,col='red')    
mtext(side=2,line=2.5,outer=F,ytext,cex=1,las=0)

mtext(side=1,line=2,outer=F,"Year",cex=1.2)

dev.off()



