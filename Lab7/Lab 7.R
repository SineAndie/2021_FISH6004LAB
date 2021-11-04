# lab 7
# install.packages("nlstools")
library("nlstools")
library("scam")
library("nlstools")
library("tidyverse");
library("readxl");
library("broom");
library("ggpmisc")

# Part 1: VonB 
k = 0.1
Linf=120
Lt=seq(0,120,length=100)

dlt = k*(Linf-Lt)

plot(Lt,dlt,xlab='L(t)',ylab='dL(t)/dt',type='l',lwd=2)

## Haddon Example
n=11
age = c(c(1,2),seq(3.3,11.3,by=1))
len =c(15.4,26.93,42.23,44.59,47.63,49.67,50.87,52.3,54.77,56.43,55.88)
print(age)
print(len)
gdata=data.frame(age=age,len=len)

VonB.fit <- nls(len ~ Linf*(1 - exp(-k*(age-ao))),
                algorithm="port",
                start = list(Linf=50, k = 0.1, ao = 0), 
                lower=c(0,0,-10),
                data=gdata)

summary(VonB.fit)
confint(VonB.fit)

VonB.predict = predict(VonB.fit)
resid = gdata$len - VonB.predict
se.resid = sqrt(sum(resid**2)/(n - 3))
resid.s = resid/se.resid

par(mfrow=c(2,1),mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(gdata$age,gdata$len,xlab='age',ylab='Length',
     type='p',lwd=2)
lines(gdata$age,VonB.predict,lwd=2,col='red')
plot(gdata$age,resid.s,xlab='age',ylab='Std. residual',
     type='p',lwd=2)
abline(h=0)

nr <- nlsResiduals(VonB.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

pdata = data.frame(age = seq(0,15,by=0.1)) 
VonB.predict = predict(VonB.fit,newdata=pdata)

ylim = range(VonB.predict)
xlim = range(pdata$age)

#bootstrapped CIs
boo <- nlsBoot(VonB.fit)
plot(boo)

par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(2,1,0))
hist(boo$coef[,1],xlab="Linf",main='')
hist(boo$coef[,2],xlab="k",main='')
hist(boo$coef[,3],xlab="ao",main='')

summary(boo)

###########  Shelton and Mangel parametriztion ############
par(mar=c(3.5,3,3,2),mgp=c(1,1,0))
VonB.fit1 <- nls(len ~ lo*exp(-k*age) + q*(1 - exp(-k*age))/k,
                 algorithm="port",
                 start = list(q=16, k = 0.1, lo = 0), 
                 lower=c(0,0,-10),upper=c(100,100,100),
                 data=gdata)
summary(VonB.fit1)

boo1 <- nlsBoot(VonB.fit1)
plot(boo1)

par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(2,1,0))
hist(boo1$coef[,1],xlab="Linf",main='')
hist(boo1$coef[,2],xlab="k",main='')
hist(boo1$coef[,3],xlab="ao",main='')

## scam fit
sfit <- scam(len ~ s(age,bs="mpi"),
             family=gaussian(link="identity"),
             data=gdata)

pdata = data.frame(age = seq(0,20,by=0.1)) 
sfit.pred = predict(sfit,newdata=pdata)       

sfit1 <- scam(len ~ s(age,bs="micv"),
              family=gaussian(link="identity"),
              data=gdata)

sfit1.pred = predict(sfit1,newdata=pdata) 

vonb.fit = predict(VonB.fit,newdata=pdata)

ylim = range(sfit.pred,sfit1.pred,vonb.fit)
xlim = range(pdata$age)

par(mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(gdata$age,gdata$len,xlab='age',ylab='Length',
     type='p',lwd=2,ylim=ylim,xlim=xlim) 
lines(pdata$age,vonb.fit,lwd=2,col='blue') 
lines(pdata$age,sfit.pred,lwd=2,col='red') 
lines(pdata$age,sfit1.pred,lwd=2,col='green')
abline(h=0,lty=2)

legend("bottom",col=c("blue","red","green"),lty=1,legend=c("VonB","Monotone increasing","Monotone increasing and concave"),bty='n',lwd=2)

#--------------------------------
# Part 2: W@A
t = 0:50
Lt = Linf*(1-exp(-k*t))
a=1;b=3
Winf = 3
a = Winf/(Linf**b)
Wt = a*(Lt)**b

par(las=1,cex.axis=1.5,cex.lab=1.5,mar=c(4,4.5,0.5,4.5))
plot(t,Lt,xlab='time, t',ylab='',type='l',lwd=2,col='red')
mtext("L(t)",side=2,line=3,las=0,cex=1.5)
par(new=TRUE)
plot(t,Wt,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",lwd=2)
axis(4)
mtext("W(t)",side=4,line=3,las=0,cex=1.5)
legend("bottomright",col=c("red","blue"),lty=1,legend=c("L(t)","W(t)"),bty='n',lwd=2)

## American plaice example
#DATA_SECTION
meanL <- read_excel("meanLWM.xlsx","meanL",na=".")%>%
  unite(index,year,div,sex,age,remove=F)%>%subset(estimate>=9.5)

meanW <- read_excel("meanLWM.xlsx","meanW",na=".")%>%
  unite(index,year,div,sex,age)%>%mutate(estimate=estimate*1000)

meanLW <- full_join(meanL,meanW,by="index")%>%na.omit()%>%
  select(div,year,age,sex,otolith,meanL="estimate.x",meanW="estimate.y",index)%>%
  mutate(cohort=year-age,
         log_weight=log(meanW),
         log_length=log(meanL))

gdata = subset(meanLW,data=meanLW,
               (sex=='female')&(div=='2J')&(year==2010))
n = length(gdata$meanW)
gdata$meanW = gdata$meanW/1000

par(mar=c(3.5,4,0.5,1),mgp=c(2,1,0))
plot(gdata$age,gdata$meanW,cex=sqrt(gdata$otolith),xlab='Age',ylab='',las=1)
mtext(side=2,line=3,"Mean Weight (Kg)")
points(gdata$age,gdata$meanW,pch=3,cex=1/sqrt(2))
ltext = paste('No. Otoliths =',min(gdata$otolith),"to",max(gdata$otolith))
legend("bottomright",bty='n',legend = ltext)


VonB.fit <- nls(meanW ~ Winf*((1 - exp(-k*(age-ao)))**3.0),
                algorithm="port",
                start = list(Winf=0.85, k = 0.2,ao=0), 
                lower=c(0.5,0.1,-5),upper=c(2,0.4,5),
                data=gdata)

summary(VonB.fit)

VonB.fit.wt <- nls(meanW ~ Winf*((1 - exp(-k*(age-ao)))**3.0),
                   algorithm="port",
                   start = list(Winf=0.7, k = 0.4,ao=2), 
                   lower=c(0,0,-5),upper=c(2,1,5),
                   data=gdata,
                   weights=otolith)

summary(VonB.fit.wt)
print(confint(VonB.fit.wt),digits=3)
confint.default(VonB.fit.wt)

# Diagnostics
pdata = data.frame(age = seq(0,20,by=0.1)) 
VonB.predict.wt = predict(VonB.fit.wt,newdata=pdata)  
VonB.predict = predict(VonB.fit,newdata=pdata)
ylim = range(VonB.predict.wt,VonB.predict,gdata$meanW)
ylim[1]=-0.5
xlim = range(pdata$age)

par(mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(gdata$age,gdata$meanW,cex=sqrt(gdata$otolith),xlab='Age',ylab='',ylim=ylim,xlim=xlim) 
mtext(side=2,line=3,"Mean Weight (Kg)")   
points(gdata$age,gdata$meanW,pch=3,cex=1/sqrt(2))
lines(pdata$age,VonB.predict.wt,lwd=2,col='red') 
lines(pdata$age,VonB.predict,lwd=2,col='blue')
legend("bottomright",col=c("blue","red"),lty=1,legend=c("Un-weighted","Weighted"),bty='n',lwd=2)
abline(h=0)

nr <- nlsResiduals(VonB.fit)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

nr <- nlsResiduals(VonB.fit.wt)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

##################
#LN 
gdata$log_meanW = log(gdata$meanW)
VonB.LNfit.wt <- nls( log_meanW~ log_Winf + 3*log((1 - exp(-k*(age-ao)))),
                     algorithm="port",
                     start = list(log_Winf=-0.4, k = 0.3,ao=0), 
                     lower=c(-2,-2,-2),upper=c(log(20),1,1.999),
                     data=gdata,
                     weights=otolith) 

summary(VonB.LNfit.wt)

pdata = data.frame(age = seq(2,14,by=0.1)) 
VonB.predict = exp(predict(VonB.LNfit.wt,newdata=pdata))

ylim = range(VonB.predict,gdata$meanW)
xlim = range(pdata$age)

par(mfrow=c(2,1),mar=c(0,4,0.5,1),oma=c(3,0,0,0),mgp=c(2,1,0))
plot(gdata$age,gdata$meanW,cex=sqrt(gdata$otolith),xlab='Age',ylab='',las=1,xaxt='n',ylim=ylim)
mtext(side=2,line=3,"Mean Weight (Kg)")
points(gdata$age,gdata$meanW,pch=3,cex=1/sqrt(2))
lines(pdata$age,VonB.predict,lwd=2,col='red')
plot(gdata$age,nr$resi2[,2],xlab='Age',ylab='',type='p',lwd=2,las=1)       
mtext(side=2,line=3,"Std. Residual")         
mtext(side=1,line=2,"Age")
abline(h=0)

nr <- nlsResiduals(VonB.LNfit.wt)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

##########
boo <- nlsBoot(VonB.LNfit.wt)
plot(boo)

## Gompertz fit
GOMP.fit.wt <- nls(meanW ~ Winf*exp(rho*exp(-k*age)),
                   algorithm="port",
                   start = list(Winf=0.7, k = 0.4,rho=-10), lower=c(0,0,-Inf),upper=c(2,1,-2),
                   data=gdata,
                   weights=otolith)
summary(GOMP.fit.wt)

boo <- nlsBoot(GOMP.fit.wt)
plot(boo)

nr <- nlsResiduals(GOMP.fit.wt)
par(mar=c(3.5,3,3,2),mgp=c(2,1,0))
plot(nr, which = 0)

pdata = data.frame(age = seq(0,20,by=0.1)) 
VonB.predict = predict(VonB.fit.wt,newdata=pdata)  
GOMP.predict = predict(GOMP.fit.wt,newdata=pdata)
ylim = range(GOMP.predict,VonB.predict,gdata$meanW)
ylim[1]=-0.5
xlim = range(pdata$age)

par(mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(gdata$age,gdata$meanW,cex=sqrt(gdata$otolith),xlab='Age',ylab='',ylim=ylim,xlim=xlim) 
mtext(side=2,line=3,"Mean Weight (Kg)")   
points(gdata$age,gdata$meanW,pch=3,cex=1/sqrt(2))
lines(pdata$age,VonB.predict,lwd=2,col='red') 
lines(pdata$age,GOMP.predict,lwd=2,col='blue')
legend("bottomright",col=c("blue","red"),lty=1,legend=c("Gompertz","VonB"),bty='n',lwd=2)
abline(h=0)

#---------------------------------
# Part 3 Multiple populations
Linf1 = 120
k1 = 0.4
Linf2 = 100
k2 = 0.3
n=50

age=0:14
pop.len1 = Linf1*(1 - exp(-k1*age))
pop.len2 = Linf2*(1 - exp(-k2*age))
plot(age,pop.len1,xlab='age',ylab='Length',
     type='l',lwd=2,las=1,col='blue')
lines(age,pop.len2,lty=1,lwd=2,col='red')
legend("bottomright",col=c("red","blue"),lty=1,legend=c("Pop2","Pop1"),bty='n',lwd=2)

rage1 = sort(rgamma(n,shape=4/2,scale=2))
rage2 = sort(rgamma(n,shape=4/2,scale=2))
Epop.len1 = Linf1*(1 - exp(-k1*rage1))
Epop.len2 = Linf2*(1 - exp(-k2*rage2))
len1 = rnorm(n,mean=Epop.len1,sd=8)
len2 = rnorm(n,mean=Epop.len2,sd=8)

ylim=c(min(len1,len2),max(len1,len2))
par(mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(age,pop.len1,xlab='age',ylab='Length',
     type='l',lwd=2,las=1,col='blue',ylim=ylim)
lines(age,pop.len2,lty=1,lwd=2,col='red')
legend("bottomright",col=c("red","blue"),lty=1,legend=c("Pop2","Pop1"),bty='n',lwd=2)
points(rage1,len1,col='blue')
points(rage2,len2,col='red')

pop = rep(1,2*n)
pop[(n+1):(2*n)]=2
age=c(rage1,rage2)
len=c(len1,len2)
gdata=data.frame(age=age,len=len,pop=pop)
gdata$ind=0
gdata$ind[gdata$pop==2]=1

VonB.fit <- nls(len ~ (Linf1 + dLinf2*ind)*(1-exp(-(k1 + dk2*ind)*age)),
                algorithm="port",data=gdata,
                start = list(Linf1=50,dLinf2=0,k1=0.2,dk2=0), 
                lower=c(0,-60,0,-0.3))
summary(VonB.fit)
confint(VonB.fit)

#### 5 pops #######;
Linf = c(80,90,100,110,120)
k = 0.3
n=50

for(i in 1:5){
  
  rage = sort(rgamma(n,shape=4/2,scale=2))
  Epop.len = Linf[i]*(1 - exp(-k*rage))
  len = rnorm(n,mean=Epop.len,sd=8)
  dati = data.frame(age=rage,length=len,pop.length = Epop.len)
  dati$pop=i
  if(i==1){dat=dati}
  if(i>1){dat=rbind(dat,dati)}
}

par(mar=c(3.5,3,0.5,1),mgp=c(2,1,0))
plot(dat$age,dat$length,xlab='age',ylab='Length',
     type='p',col=dat$pop)
for(i in 1:5){
  ind = dat$pop==i 
  lines(dat$age[ind],dat$pop.length[ind],lty=1,lwd=2,col=i)
}

gfit = function(map,length,age,Linf_main,k_main,Linf_dev,k_dev){
  Linf_dev = c(0,Linf_dev) 
  k_dev = c(0,k_dev)
  Linf = Linf_main + Linf_dev 
  k = k_main + k_dev
  pred = Linf[map]*(1 - exp(-k[map]*age)) 
  res = length - pred;
  return(res)}

n.pop = length(unique(dat$pop))
dat$map = as.numeric(as.factor(dat$pop))

start.parm = list(Linf_main=100,k_main=0.3,
                  Linf_dev=rep(0,n.pop-1),k_dev=rep(0,n.pop-1))

fit <- nls( ~ gfit(map,length,age,Linf_main,k_main,Linf_dev,k_dev),
            data=dat,start = start.parm,  algorithm="port",control=list(maxiter=500),
            lower=c(30,0.1,rep(-Inf,n.pop-1),rep(-0.3,n.pop-1)),
            upper=c(200,0.5,rep(Inf,n.pop-1),rep(0.3,n.pop-1)))

summary(fit)
confint(fit)   

