# Lab 9 S-R relationship
## 4VWX cod Example
fname = c('SR.dat')
sr.data = read.table(fname,header=T)

plot(sr.data$year,sr.data$ssb,xlab='Year',ylab='Parental stock size, S',
     type='l')

plot(sr.data$year,sr.data$rec,xlab='Year',ylab='Recruitment, R',
     type='l')

plot(sr.data$ssb,sr.data$rec,xlab='Parental stock size, S',
     ylab='Recruitment, R',type='p',xlim=c(0,max(sr.data$ssb)))


sr.data$logrec = log(sr.data$rec)
sr.data$logssb = log(sr.data$ssb)

#choose realistic starting values
init.Rmax = 150
init.S50 = 50

init.alpha= init.Rmax
init.beta = init.S50

BH.fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb),
              algorithm="port",lower=c(0,0),
              data=sr.data,
              start = list(beta = init.beta,alpha = init.alpha))

summary(BH.fit)

BH.profile = confint(BH.fit)

n = length(sr.data$ssb)
BH.predict = predict(BH.fit)
resid = sr.data$logrec - BH.predict
se.resid = sqrt(sum(resid**2)/(n - 2))
resid.s = resid/se.resid

ssb.pred = 0:150
BH.predict = predict(BH.fit,list(ssb=ssb.pred,logssb=log(ssb.pred)))
parm.est = coef(BH.fit)

plot(sr.data$ssb,sr.data$rec,xlab='Parental stock size, S',
     ylab='Recruitment, R',type='p',xlim=c(0,max(sr.data$ssb)))

lines(ssb.pred,exp(BH.predict),lwd=2,col='green')
p1 = min(sr.data$rec)
points(parm.est[1],p1,pch=18,cex=2,col='red')
points(0.1,parm.est[2],pch=18,cex=2,col='blue')
arrows(BH.profile[1,1],p1,BH.profile[1,2],p1,col='red',lwd=2,code=3)
arrows(0.1,BH.profile[2,1],0.1,BH.profile[2,2],col='blue',lwd=2,code=3)


par(mfrow=c(2,1),mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(sr.data$ssb,sr.data$rec,xlab='SSB',ylab='Recruitment',
     type='p',lwd=2)
lines(ssb.pred,exp(BH.predict),lwd=2,col='red')
lines(ssb.pred,exp(BH.predict + (se.resid**2)/2),lwd=2,col='purple')

parm.s = log(rev(coef(BH.fit)))
y = sr.data$rec
x = sr.data$ssb

legend('topleft',lty=1,lwd=2,col=c('red','purple'),
       c("LN biased","LN corrected"),bty='n',cex=0.8,inset=-0.05)

plot(sr.data$ssb,resid.s,xlab='SSB',ylab='Std. residual',
     type='p',lwd=2)
abline(h=0)

#install.packages("nlstools")
library("nlstools")
nr <- nlsResiduals(BH.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

par(mfrow=c(2,1),mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(sr.data$year,resid.s,xlab='Year',ylab='Std. residual',type='b')
abline(h=0,lty=2)

acf(resid.s,plot=TRUE)

test.nlsResiduals(nr)

boo <- nlsBoot(BH.fit)
plot(boo)

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(2,1,0))
hist(boo$coef[,1],xlab="beta, S50",main='')
hist(boo$coef[,2],xlab="alpha, Rmax",main='')

summary(boo)

################## Ricker ###########################
init.beta = 1/(4.44*init.S50)
init.alpha= init.Rmax*init.beta*exp(1)

RK.fit <- nls(logrec ~ log(alpha) + logssb - beta*ssb,
              algorithm="port",lower=c(0,0),data=sr.data,
              start = list(beta = init.beta,alpha = init.alpha))

summary(RK.fit)

RK.profile = confint(RK.fit)

n = length(sr.data$ssb)
RK.predict = predict(RK.fit)
resid = sr.data$logrec - RK.predict
se.resid = sqrt(sum(resid**2)/(n - 2))
resid.s = resid/se.resid

RK.predict = predict(RK.fit,list(ssb=ssb.pred,logssb=log(ssb.pred)))

plot(sr.data$ssb,sr.data$rec,xlab='Parental stock size, S',
     ylab='Recruitment, R',type='p',xlim=c(0,max(sr.data$ssb)))
lines(ssb.pred,exp(BH.predict),lwd=2,col='green')
lines(ssb.pred,exp(RK.predict),lwd=2,col='red')
legend('topleft',lty=1,lwd=2,col=c('red','green'),
       c("RK","BH"),bty='n',cex=0.8)


par(mfrow=c(2,1),mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(sr.data$ssb,sr.data$rec,xlab='SSB',ylab='Recruitment',
     type='p',lwd=2)
lines(ssb.pred,exp(RK.predict),lwd=2,col='red')
lines(ssb.pred,exp(RK.predict + (se.resid**2)/2),lwd=2,col='purple')

parm.s = log(rev(coef(RK.fit)))
y = sr.data$rec
x = sr.data$ssb

legend('topleft',lty=1,lwd=2,col=c('red','purple'),
       c("LN biased","LN corrected"),bty='n',cex=0.8,inset=-0.05)

plot(sr.data$ssb,resid.s,xlab='SSB',ylab='Std. residual',
     type='p',lwd=2)
abline(h=0)

nr <- nlsResiduals(RK.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

par(mar=c(3.5,3,1,1),mgp=c(2,1,0))
plot(sr.data$year,resid.s,xlab='Year',ylab='Std. residual',type='b')
abline(h=0,lty=2)

test.nlsResiduals(nr)

boo <- nlsBoot(RK.fit)
plot(boo)

boo.Rmax = boo$coef[,2]/(boo$coef[,1]*exp(1));

ricker <-function(parm,s){
  alpha <- parm[2]
  beta <- parm[1];
  ret = alpha*s*exp(-beta*s);
  return(ret)
}

fS50 <- function(parm){
  hRmax = 0.5*parm[2]/(parm[1]*exp(1));
  temp = uniroot(function(x) hRmax-ricker(parm,x),
                 lower = 0, upper = 1/parm[1])
  return(temp$root)
}

boo.S50 = apply(boo$coef,1,fS50);

new.boo = boo
new.boo$coef =cbind(boo$coef,boo.Rmax,boo.S50)
colnames(new.boo$coef) = c(colnames(boo$coef),"Rmax","S50")
qRmax = quantile(boo.Rmax,probs=c(0.5,0.025,0.975))
qS50 = quantile(boo.S50,probs=c(0.5,0.025,0.975))
new.boo$bootCI = rbind(boo$bootCI,qRmax,qS50)
rownames(new.boo$bootCI) = c(rownames(boo$bootCI),"Rmax","S50")

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(2,1,0))
hist(new.boo$coef[,4],xlab="S50",main='')
hist(new.boo$coef[,3],xlab="Rmax",main='')

options(digits=3)
summary(new.boo)

###  hockey stick  ############################
init.alpha= 0.5*init.Rmax/init.S50
init.delta = init.S50
gam2by4 = 0.1/4
delta.max = max(x) 
delta.min = min(x)

HS.fit <- nls(logrec ~ log(alpha) + log(ssb + sqrt(delta**2 + gam2by4) - sqrt((ssb-delta)**2 + gam2by4)),
              data=sr.data,
              start = list(alpha = init.alpha,delta = init.delta),
              algorithm="port",lower=c(0,delta.min),upper=c(25,delta.max))
summary(HS.fit)

gam2by4 = 500/4

HS.fit1 <- nls(logrec ~ log(alpha) + log(ssb + sqrt(delta**2 + gam2by4) - sqrt((ssb-delta)**2 + gam2by4)),
               data=sr.data,
               start = list(alpha = init.alpha,delta = init.delta),
               algorithm="port",lower=c(0,delta.min),upper=c(25,delta.max))

HS.parm = coef(HS.fit)  
HS.profile = confint(HS.fit)

n = length(sr.data$ssb)
HS.predict = predict(HS.fit)  
resid = sr.data$logrec - HS.predict
se.resid = sqrt(sum(resid**2)/(n - 2))
resid.s = resid/se.resid

HS.predict = predict(HS.fit,list(ssb=ssb.pred,logssb=log(ssb.pred))) 
HS.predict1 = predict(HS.fit1,list(ssb=ssb.pred,logssb=log(ssb.pred)))

par(mar=c(3,3,0.5,1),mgp=c(2,1,0))
plot(sr.data$ssb,sr.data$rec,xlab='Parental stock size, S',
     ylab='Recruitment, R',type='p',xlim=c(0,max(sr.data$ssb)))

lines(ssb.pred,exp(BH.predict),lwd=2,col='green')
lines(ssb.pred,exp(RK.predict),lwd=2,col='red')   
lines(ssb.pred,exp(HS.predict),lwd=2,col='blue')  
lines(ssb.pred,exp(HS.predict1),lwd=2,lty=2,col='blue')
legend('topleft',lwd=2,col=c('red','green','blue','blue'),lty=c(1,1,1,2),
       c("RK","BH","HS","HS gam=500"),bty='n',cex=0.8)


HS.boo <- nlsBoot(HS.fit)

par(mar=c(2,2,0.5,0.5),mgp=c(3,1,0))
plot(HS.boo)

########### SCAM #####################
require("scam")
sr.data$offset = log(sr.data$ssb)

scam.fit <- scam(logrec ~ s(ssb,k=5,bs="mpd",m=2) + offset(offset),
                 family=gaussian(link="identity"),data=sr.data,sp=0.000001)

summary(scam.fit)

scam.predict = predict(scam.fit,list(ssb=ssb.pred,offset=log(ssb.pred)))         

par(mfrow=c(2,1),oma=c(2,0,0,0),mar=c(0,3,0.5,1),mgp=c(2,1,0))

plot(sr.data$ssb,sr.data$rec,xlab='S',
     ylab='R',type='p',xlim=c(0,max(sr.data$ssb)),xaxt='n')
lines(ssb.pred,exp(BH.predict),lwd=2,col='green')
lines(ssb.pred,exp(RK.predict),lwd=2,col='red')   
lines(ssb.pred,exp(HS.predict),lwd=2,col='blue')  
lines(ssb.pred,exp(scam.predict),lwd=2,col='purple')
legend('topleft',lwd=2,col=c('red','green','blue','purple'),lty=c(1,1,1,1),
       c("RK","BH","HS","SCAM"),bty='n',cex=0.8)

plot(sr.data$ssb,sr.data$rec/sr.data$ssb,xlab='S',
     ylab='R/S',type='p',xlim=c(0,max(sr.data$ssb)))
lines(ssb.pred,exp(BH.predict)/ssb.pred,lwd=2,col='green')
lines(ssb.pred,exp(RK.predict)/ssb.pred,lwd=2,col='red')   
lines(ssb.pred,exp(HS.predict)/ssb.pred,lwd=2,col='blue')  
lines(ssb.pred,exp(scam.predict)/ssb.pred,lwd=2,col='purple')


######## BH scam #################
#first
sr.data$ssb_inv = 1/sr.data$ssb

BH.fit <- glm(rec ~ ssb_inv, family = Gamma(link="inverse"), data = sr.data)  
summary(BH.fit) 
BH.predict = predict(BH.fit,list(ssb_inv=1/ssb.pred),type = "response")

scamBH.fit <- scam(rec ~ s(ssb_inv,k=5,bs="mpi",m=2),
                   family=Gamma(link="inverse"),data=sr.data)

summary(scamBH.fit)

scamBH.predict = predict(scamBH.fit,list(ssb_inv=1/ssb.pred),type = "response")         

par(mar=c(3,3,0.5,1),mgp=c(2,1,0))

plot(sr.data$ssb,sr.data$rec,xlab='S',
     ylab='R',type='p',xlim=c(0,max(sr.data$ssb)),xaxt='n')
lines(ssb.pred,BH.predict+1,lwd=2,col='green')
lines(ssb.pred,exp(RK.predict)+1,lwd=2,col='red')   
lines(ssb.pred,exp(scam.predict),lwd=2,col='blue') 
lines(ssb.pred,scamBH.predict,lwd=2,col='black')
legend('topleft',lwd=2,col=c('red','green','blue','black'),lty=c(1,1,1,1),
       c("RK (jittered)","Gamma BH  (jittered)","SCAM_RK","SCAM_BH"),bty='n',cex=0.8)

##remember correlated errors 
#### BH with AR1 errors
BH.fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb),
              algorithm="port",lower=c(0,0),data=sr.data,
              start = list(beta = init.beta,alpha = init.alpha))


BH.arfit <- gnls(logrec ~ log(alpha) + logssb - log(beta + ssb),data=sr.data,
                 start = list(beta = init.beta,alpha = init.alpha), 
                 correlation = corAR1(form=~year))

summary(BH.arfit) 

BH.predict = predict(BH.fit,list(ssb=ssb.pred,logssb=log(ssb.pred)))
BH.arpredict = predict(BH.arfit,list(ssb=ssb.pred,logssb=log(ssb.pred)))

resid.s=residuals(BH.arfit, type='normalized')


par(mfrow=c(2,1),oma=c(3,0,0,0),mar=c(0,3,0.5,1),mgp=c(2,1,0))

plot(sr.data$ssb,sr.data$rec,xlab='SSB',ylab='Recruitment',
     type='p',lwd=2,xaxt='n')
lines(ssb.pred,exp(BH.predict),lwd=2,col='green')
lines(ssb.pred,exp(BH.arpredict),lwd=2,col='red')

legend('topleft',lwd=2,col=c('green','red'),lty=c(1,1),
       c("BH","BH AR"),bty='n',cex=0.8)

plot(sr.data$ssb,resid.s,xlab='SSB',ylab='Std. residual',type='p',lwd=2)
mtext(side=1,line=2,outer=T,'SSB')
abline(h=0)


par(mfrow=c(2,1),mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(sr.data$year,resid.s,xlab='Year',ylab='Std. residual',type='b')
abline(h=0,lty=2)

acf(resid.s,plot=TRUE)

par(mar=c(3,3,1,1),mgp=c(2,1,0))
qqnorm(resid.s) 
qqline(resid.s)

shapiro.test(resid.s)










