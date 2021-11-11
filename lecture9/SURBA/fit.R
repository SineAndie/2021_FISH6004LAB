setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week8\\F3LNO_SURBA")

source("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week2\\spay.txt")
 
require(lattice)
library(latticeExtra)
require(reshape2)  
library(lme4)          
library(cAIC4)  
library(xtable)  
library(ggplot2)  
library(psych)    
library(stargazer)  
library(numDeriv)

load('F3LNO.RData')
head(vdat,n=10) 
head(mdat,n=10)

vdat$YC = vdat$Year-vdat$Age
vdat$log.index = NA
ind = vdat$wt==1
vdat$log.index[ind] = log(vdat$index[ind])

my.padding <- list(layout.heights = list( 
                        top.padding = 0, 
                        main.key.padding = 0, 
                        key.axis.padding = 0, 
                        axis.xlab.padding = 0, 
                        xlab.key.padding = 0, 
                        key.sub.padding = 0), 
                layout.widths = list( 
                        left.padding = 1,
                        key.ylab.padding = 0, 
                        ylab.axis.padding = 1, 
                        axis.key.padding = 0, 
                        right.padding = 0) 
                ) 

colv = colors()[seq(26,by=5,length=13)]

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Index"), cex=1)
stripttl <- list(cex=0.75,font=1) 
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
my.key <- list(x = .96, y = .96, corner = c(1, 1),
            size = 2,between=0.5,
            text = list(as.character(0:12),cex=0.6),
            lines = list(type = c("l"), col = colv,
                         lwd = 2, lty = 1,span=0.25))
 
P1 = xyplot(index~Year,groups=Age,cex=1,ylim=c(0,6),
data=vdat, type=c("l"),lty=1,xlab=xttl,lwd=2,
ylab=yttl, strip.text=stripttl, las=1,
key=my.key,col=colv,par.settings = my.padding,
par.strip.text = list(cex=0.6),
) 
print(P1)

gname=c('RV.jpeg')
jpeg(file=gname,width=4,height=4,units='in',res=300)
par(mar=c(2.5,3,0,0),mgp=c(2,0.75,0))
print(P1) 
dev.off()



vdat$log.index.dev=NA
vdat$log.index.dev[ind]=resid(glm(log.index ~ factor(Age) -1,data=vdat))

my.key <- list(corner = c(1,0),
            size = 2,between=0.5,
            text = list(as.character(0:12),cex=0.6),
            lines = list(type = c("l"), col = colv,
                         lwd = 2, lty = 1,span=0.25))

xttl <- list(c("Cohort"), cex=1)
yttl <- list(c("Log Index Deviation"), cex=1)
P2 = xyplot(log.index.dev~YC,groups=Age,cex=1,
data=vdat, type=c("l"),lty=1,xlab=xttl,lwd=2,
ylab=yttl, strip.text=stripttl, las=1,
key=my.key,col=colv,par.settings = my.padding,
par.strip.text = list(cex=0.6),
) 
print(P2)

gname=c('log_RV_dev.jpeg')
jpeg(file=gname,width=4,height=4,units='in',res=300)
par(mar=c(2.5,3,0,0),mgp=c(2,0.75,0))
print(P2) 
dev.off() 

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Age"), cex=1)
P3 = xyplot(Age~Year,cex=1,
data=vdat, type=c("p"),xlab=xttl,pch=16*vdat$wt,
ylab=yttl, strip.text=stripttl, las=1,
par.settings = my.padding,
par.strip.text = list(cex=0.6),
) 
print(P3)

gname=c('RV_wt.jpeg')
jpeg(file=gname,width=4,height=4,units='in',res=300)
par(mar=c(2.5,3,0,0),mgp=c(2,0.75,0))
print(P3) 
dev.off() 



scale=1
 
sp.dat <- t(spay(t(mdat[,2:14]))) 
syear <- rep(mdat[,1],length(0:12))
sage <- rep(0:12,each=length(mdat[,1]))
vspay <- as.vector(unlist(sp.dat))

cohort.ps=c(0.02,0.01,0.03,0.01)
y.adj = 0.2
ylim=range(sage);
ylim[2]=ylim[2]+y.adj
xlim=range(syear);
xlim[2]=xlim[2]+0.2

gname <- 'spay.jpeg'
jpeg(file=gname,width=4,height=4,units='in',res=300)

par(mar=c(3,3,1.2,0.2),cex.axis=0.7)

ind = !is.na(vspay)
bp(syear[ind],sage[ind],vspay[ind],cohort.ps=cohort.ps, ylim=ylim, xlim=xlim, 
     las=1, xlab='', ylab='', main="", scale=scale)
points(syear[ind],sage[ind],cex=sqrt(abs(vspay[ind]))*scale, col='black', pch=1)
mtext(side=3,outer=F,line=.1,c('3NO Cod Fall RV (ages 0-12) SPAY'),cex=0.8) 
mtext(side=2,outer=F,line=2,'Age')  
mtext(side=1,outer=F,line=2,'Year')
  
dev.off()


## create some useful objects

uage = unique(vdat$Age)
uyear = min(vdat$Year):max(vdat$Year)
A = length(uage)
Y = length(uyear)
## the vectorized model dimensions to create an index observation map
pdat = data.frame(
  Age = as.vector(matrix(uage,nrow=Y,ncol=A,byrow=T)),
  Year = as.vector(matrix(uyear,nrow=Y,ncol=A,byrow=F))
)
pdat$imap = 1:nrow(pdat)

vdat.nz = subset(vdat,wt==1); #use data with wts=1
temp1 = paste(pdat$Year,"_",pdat$Age,sep='')
temp2 = paste(vdat.nz$Year,"_",vdat.nz$Age,sep='')
vdat.nz$imap = pdat$imap[temp1 %in% temp2]



persp(mdat[,1],uage,as.matrix(mdat[,2:14]), 
      main="Survey Numbers-at-age Surface",
      xlab='Year',ylab='Age',zlab='Number',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
      
      
## First trial values for survey index q_age
d_o = exp(-0.6*(0:(A-1)))
mean_age = aggregate(vdat$index,list(age=vdat$Age),mean)$x
q_age = mean_age/d_o
q_age=q_age/max(q_age)
plot(uage,q_age)

q_age = c(0.005,0.3,1,1,1,1,1,1,1,1,1,1,1) 
Qm = matrix(q_age,nrow=Y,ncol=A,byrow=T)


gname <- 'q.jpeg'
jpeg(file=gname,width=2,height=2,units='in',res=300)
  par(mar=c(3,3,0.2,0.2))
  plot(uage,q_age,type='l',lwd=2,xlab='',ylab='')
  mtext(side=2,outer=F,line=2,'Index Q')  
  mtext(side=1,outer=F,line=2,'Age') 
dev.off()

Npop = function(Nfirst,No,Z){
  Nmatrix = matrix(NA,nrow=Y,ncol=A)
  next.age = 2:A
  Nmatrix[,1] = No
  Nmatrix[1,next.age] = Nfirst
  for(y in 2:Y){
    Nmatrix[y,next.age] = Nmatrix[y-1,next.age-1]*exp(-Z[y-1,next.age-1]);
  }
  return(Nmatrix)
}

      
pred_index = function(logNo,logNfirst,sparm,logf){
  No = exp(logNo)
  Nfirst=exp(logNfirst)
  f = exp(logf)
  age_diff = uage-6
  s_age = exp(sparm[1]*age_diff + sparm[2]*age_diff*age_diff)  
  Z = f %o% s_age
  N = Npop(Nfirst,No,Z)
  B = wtm*N
  SSB = mat*B
  logR = log(Qm*N*exp(-0.75*Z))
  vN = as.vector(N)   
  vB = as.vector(B) 
  vSSB = as.vector(SSB) 
  vZ = as.vector(Z) 
  vlogR = as.vector(logR)
  return(cbind(vlogR,vN,vB,vSSB,vZ))
}

N.start = mdat[,2:14]/Qm
Nfirst = apply(N.start,2,mean)[2:A]
No = rep(mean(N.start[,1]),Y)
sparm=c(0,0) 
age_diff = uage-6
s_age = exp(sparm[1]*age_diff + sparm[2]*age_diff*age_diff)
f_year = rep(0.4,Y)
Z = f_year %o% s_age
Ninit = Npop(Nfirst,No,Z)


start.parms = list(
  logNo = log(No),
  logNfirst = log(Nfirst),
  sparm=c(0,0),
  logf = log(f_year)
)

test = function(x){  
  logNo = x[1:27];
  logNfirst = x[28:39];
  sparm = x[40:41];
  logf = x[42:length(x)];
  return(pred_index(logNo,logNfirst,sparm,logf)[vdat.nz$imap,1])
 }

x = unlist(start.parms) 
test.jac = jacobian(test, x)
std.test.jac = sqrt(apply(test.jac,2,var))
plot(std.test.jac,type='l',ylim=c(0,0.2))  
abline(h=0,col='red')
x[which(std.test.jac<1e-6)]

plot(test(x),vdat.nz$log.index,xlab='init pred',ylab='observed')
abline(a=0,b=1,col='red')

fit = function(x){
  return(sum((vdat.nz$log.index-test(x))**2))
 }
x.start=optim(x,fit,control=list(trace=1))$par 

start.parms = list( 
  logNo = x.start[1:27],
  logNfirst = x.start[28:39],
  sparm = x.start[40:41],
  logf = x.start[42:length(x)]
)
start.parms$logf[start.parms$logf<log(0.2)]=log(0.3)
 
lower = list(
  logNo = rep(-Inf,length(start.parms$logNo)),
  logNfirst = rep(-Inf,length(start.parms$logNfirst)),
  sparm = rep(-Inf,length(start.parms$sparm)),
  logf = rep(log(0.2),length(start.parms$logf))
)

upper = list(
  logNo = rep(Inf,length(start.parms$logNo)),
  logNfirst = rep(Inf,length(start.parms$logNfirst)),
  sparm = rep(Inf,length(start.parms$sparm)),
  logf = rep(Inf,length(start.parms$logf))
)
lower=unlist(lower)
upper=unlist(upper) 

surba.fit <- nls(log.index ~ pred_index(logNo,logNfirst,sparm,logf)[imap,1], 
  algorithm="port",data=vdat.nz,start = start.parms,
  control=list(maxiter=10000),lower=lower,upper=upper)
  
sum.surba = summary(surba.fit)


stargazer(sum.surba[11], type = "html",  out="out.doc",align = TRUE,   
 title="Fall Survey 3LNO cod SURBA", single.row=TRUE)
# ci=TRUE, ci.level=0.95)

x = coef(surba.fit)
logNo = x[1:27]
logNfirst = x[28:39]
sparm = x[40:41]
logf = x[42:length(x)]
  
No = exp(logNo)
Nfirst=exp(logNfirst)
f = exp(logf)
age_diff = uage-6
s = exp(sparm[1]*age_diff + sparm[2]*age_diff*age_diff)  
Z = f %o% s
N = Npop(Nfirst,No,Z)

pred.dat = data.frame(pred_index(logNo,logNfirst,sparm,logf))
pred.dat$Year=rep(uyear,13)
pred.dat$Age = rep(0:12,each=Y)
pred.dat$F = pred.dat$vZ - 0.2

temp = subset(pred.dat,(Age>=4)&(Age<=6))
temp1 = aggregate(temp$F,list(Year=temp$Year),mean)

gname <- 'aveF.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)
  par(mar=c(3,3,0.2,0.2))
  plot(temp1$Year,temp1$x,type='l',lwd=2,xlab='',ylab='',xlim=c(1959,2016))
  mtext(side=2,outer=F,line=2,'F = Z-0.2 (Ave, Ages 4-6)')  
  mtext(side=1,outer=F,line=2,'Year') 
dev.off()

temp1 = aggregate(pred.dat$vSSB,list(Year=pred.dat$Year),sum)

assmnt.ssb = c(40088,30523,14503,7309,4905,5998,8076,8496,8598,8576,7577,7906,
7780,9646,7447,8110,7586,6831,8224,8956,9047,11594,17749,21610,23304,24349,23204)

rel.ssb = assmnt.ssb*mean(temp1$x)/mean(assmnt.ssb)

ylim=range(rel.ssb,temp1$x)

gname <- 'ssb.jpeg'
jpeg(file=gname,width=3,height=3,units='in',res=300)
  par(mar=c(3,3,0.2,0.2))
  plot(temp1$Year,temp1$x,type='l',lwd=2,xlab='',ylab='',ylim=ylim)
  lines(temp1$Year,rel.ssb,col='red',lwd=2)
  mtext(side=2,outer=F,line=2,'Relative SSB')  
  mtext(side=1,outer=F,line=2,'Year') 
  legend("topright",col=c('black','red'),lwd=2,bty='n',legend=c('SURBA','Assessment'))
dev.off()


ssb.matrix = matrix(pred.dat$vSSB,nrow=Y,ncol=A,byrow=FALSE)

gname = "SSB_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
  persp(uyear,uage,ssb.matrix, 
      main="SSB-at-age Surface",
      xlab='Year',ylab='Age',zlab='Weight',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
dev.off() 

gname = "Nold_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(uyear,uage[6:13],N[,6:13], 
      main="Numbers-at-age Surface",
      xlab='Year',ylab='Age',zlab='Number',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
dev.off()

gname = "N_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(uyear,uage,N, 
      main="Numbers-at-age Surface",
      xlab='Year',ylab='Age',zlab='Number',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
dev.off() 

gname = "Z_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(uyear,uage,Z, 
      main="Z-at-age Surface",
      xlab='Year',ylab='Age',zlab='Z',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
dev.off()

Sm=exp(-Z)
gname = "Surv_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
persp(uyear,uage,Sm, 
      main="Survival-at-age Surface",
      xlab='Year',ylab='Age',zlab='Survival',
      col="lightgreen", 
      theta=130, phi=40,
      r=50,
      d=0.1,
      expand=0.5,
      ltheta=90, lphi=180,
      shade=0.25,
      ticktype="detailed",
      nticks=5) # produces the 3-D plot
dev.off()

ylim = range(100*Sm[,c(1,6,13)])

gname = "Surv_ts.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)
  par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(100*Sm[,1]~uyear, ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Year") 
  mtext(side=2,line=2.3,"Survival (%)")
  lines(uyear,100*Sm[,6],col='blue',lwd=2)
  lines(uyear,100*Sm[,13],col='red',lwd=2)
  legend("bottom",col=c('black','blue','red'),lwd=2,legend=uage[c(1,6,13)],bty='n')
dev.off()

pcatch = wtm*N*(1-exp(-Z))*(abs(Z-0.2))/Z
rownames(pcatch)=uyear
colnames(pcatch)=uage
pland = apply(pcatch,1,sum)

assmnt.catch = c(28846,29454,12752,10646,2702,172,
174,383,547,919,1050,1310,2194,4870,934,
724,600,848,923,1083,946,867,734,1113,734,586,666)

rel.catch = assmnt.catch*mean(pland)/mean(assmnt.catch)

ylim=range(rel.catch,pland)

gname <- 'catch_trend.jpeg'
jpeg(file=gname,width=3,height=4,units='in',res=300)
#  par(mar=c(3,3,0.2,0.2))
  par(mfcol=c(2,1),oma=c(0,1,0.5,0.5),mar=c(3,2,0,0),las=1)
  plot(uyear,pland,type='l',lwd=2,xlab='',ylab='',ylim=ylim)
  lines(uyear,rel.catch,col='red',lwd=2)
  mtext(side=2,outer=F,line=2,'Relative Catch',las=0)  
  mtext(side=1,outer=F,line=2,'Year') 
  legend("topright",col=c('black','red'),lwd=2,bty='n',legend=c('SURBA','Assessment'))
  
  
  plot(log(rel.catch),log(pland),type='p',xlab='',ylab='')
  abline(a=0,b=1,col='red')
  mtext(side=2,outer=F,line=2,'SURBA log relative Catch',las=0)  
  mtext(side=1,outer=F,line=2,'Assessment log relative Catch')
  ltext = paste('corr=',round(cor(rel.catch,pland),digits=3),' log corr=',        
  round(cor(log(rel.catch),log(pland)),digits=3))
  legend('bottomright',legend=ltext,bty='n',cex=0.5) 
  
dev.off()
  
    

fit=vdat.nz
fit$survey="Fall 3LNO"
fit$resid = residuals(surba.fit,type='pearson')
fit$pred = predict(surba.fit)

gname <- "sres.jpeg"
jpeg(file=gname,width=4,height=5,units='in',res=300)
par(mfcol=c(4,1),oma=c(0,3,1,1),mar=c(3,2,0,2),las=1)

ind <- !is.na(fit$pred)
x <- fit$Year[ind]
y <- fit$resid[ind]

plot(x,y,xlab="",ylab="",type='n')
text(x,y,fit$Age[ind])
abline(h=0,lty=1)
mres <- tapply(y,x,"mean")
lines(as.numeric(names(mres)),mres,lty=1,lwd=1,col='red')
mtext(side=4,line=0.5,outer=F,'Year',las=0)

ind1 = ind
x1 = fit$YC[ind1] 
y1 <- fit$resid[ind1]
plot(x1,y1,xlab="",ylab="",type='n')   
text(x1,y1,fit$Age[ind])
abline(h=0,lty=1)  
mres <- tapply(y1,x1,"mean")
lines(sort(unique(x1)),mres,lty=1,lwd=1,col='red')
mtext(side=4,line=0.5,outer=F,'Cohort',las=0)

plot(fit$Age[ind],y,xlab="",ylab="",pch=3)
abline(h=0,lty=2)
mtext(side=4,line=0.5,outer=F,'Age',las=0)
ma = tapply(y,fit$Age[ind],mean)
lines(as.numeric(names(ma)),ma,lty=1)

plot(fit$pred[ind],y,xlab="",ylab="",pch=3)
abline(h=0,lty=1)
mtext(side=4,line=0.5,outer=F,'Expected',las=0)

i = "Fall 3LNO"
ytext <- paste("Standardized",i,"Residuals")
mtext(side=2,ytext,line=1,outer=T,las=0,cex=1.2)

dev.off()



fit$colr = 'deepskyblue'
fit$colr[fit$resid>0]='firebrick2'

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year , data=fit, 
       ylab='Age',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 3*sqrt(abs(fit$resid)/pi), fill.color = fit$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(h=1:6, v=c(2000,2005,2010,2015), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
print(resid.plot1)
                  
gname = c("resid_matrix_bubbles.jpeg")
jpeg(file=gname,width=4,height=3,units='in',res=300)
  print(resid.plot1)
dev.off()

##Projections

F = Z-0.2
ave.F = apply(F[(Y-2):Y,],2,mean)
sel_a = ave.F/max(ave.F)


vdat.all = vdat
wtm.all=wtm
mat.all=mat
Qm.all=Qm

retro.years=seq(2016,2010,by=-1)

retro = vector("list", length(retro.years))

for(i in 1:length(retro.years)){

ry=retro.years[i]
vdat = subset(vdat.all,Year<=ry)

uyear = min(vdat$Year):max(vdat$Year)
Y = length(uyear)
pdat = data.frame(
  Age = as.vector(matrix(uage,nrow=Y,ncol=A,byrow=TRUE)),
  Year = as.vector(matrix(uyear,nrow=Y,ncol=A,byrow=FALSE))
)
pdat$imap = 1:nrow(pdat)  
wtm=wtm.all[1:Y,]
mat=mat.all[1:Y,]  
Qm=Qm.all[1:Y,]

vdat.nz = subset(vdat,wt==1); #use data with wts=1
temp1 = paste(pdat$Year,"_",pdat$Age,sep='')
temp2 = paste(vdat.nz$Year,"_",vdat.nz$Age,sep='')
vdat.nz$imap = pdat$imap[temp1 %in% temp2]

start.parms = list( 
  logNo = logNo[1:Y],
  logNfirst = logNfirst,
  sparm = sparm,
  logf = logf[1:Y]
)
start.parms$logf[start.parms$logf<log(0.2)]=log(0.3)
 
lower = list(
  logNo = rep(-Inf,length(start.parms$logNo)),
  logNfirst = rep(-Inf,length(start.parms$logNfirst)),
  sparm = rep(-Inf,length(start.parms$sparm)),
  logf = rep(log(0.2),length(start.parms$logf))
)

upper = list(
  logNo = rep(Inf,length(start.parms$logNo)),
  logNfirst = rep(Inf,length(start.parms$logNfirst)),
  sparm = rep(Inf,length(start.parms$sparm)),
  logf = rep(Inf,length(start.parms$logf))
)
lower=unlist(lower)
upper=unlist(upper) 

surba.fit.retro <- nls(log.index ~ pred_index(logNo,logNfirst,sparm,logf)[imap,1], 
  algorithm="port",data=vdat.nz,start = start.parms,
  control=list(maxiter=10000),lower=lower,upper=upper)

parm.est = coef(surba.fit.retro)  
p1 = parm.est[substr(names(parm.est),1,5)=='logNo'] 
p2 = parm.est[substr(names(parm.est),1,5)=='logNf'] 
p3 = parm.est[substr(names(parm.est),1,5)=='sparm'] 
p4 = parm.est[substr(names(parm.est),1,4)=='logf']
  
pred.dat = data.frame(pred_index(p1,p2,p3,p4))
pred.dat$Year=rep(uyear,13)
pred.dat$Age = rep(0:12,each=Y)
pred.dat$F = pred.dat$vZ - 0.2

retro[[i]]=list(pred.dat=pred.dat,retro.year=ry)

}


retro.ssb=lapply(retro,function(x){aggregate(x$pred.dat$vSSB,list(Year=x$pred.dat$Year),sum)})

ylim=range(unlist(lapply(retro.ssb,function(y){y$x[y$Year>=1998]})))
xlim=c(1998,2016)

gname <- 'ssb_retro.jpeg'
jpeg(file=gname,width=3,height=3,units='in',res=300)

  par(mar=c(3,3,0.2,0.2))
  plot(xlim,ylim,type='n',lwd=2,xlab='',ylab='',ylim=ylim)
  for(i in 1:length(retro.years)){
    x = retro.ssb[[i]]$Year
    y = retro.ssb[[i]]$x
    lines(x,y,lwd=2)
    points(x[length(x)],y[length(x)],cex=2)
  }
  mtext(side=2,outer=F,line=2,'Relative SSB')  
  mtext(side=1,outer=F,line=2,'Year') 
  
dev.off()





