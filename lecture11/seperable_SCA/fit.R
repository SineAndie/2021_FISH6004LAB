setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week10\\3NOcod")

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
                        right.padding = 1) 
                ) 

require(lattice)
library(latticeExtra)
require(reshape2)  
library(lme4)    
library(ggplot2)  
library(psych)    
library(stargazer)  
library(numDeriv)   
library(splines)

dname = '..\\..\\week9\\3NOcod_ADAPT\\data\\'
load(paste(dname,'3LNO.RData',sep=''))

# get log indices and wt=0 for index=0
indices$YC = indices$Year-indices$Age 
#set fraction of year survey takes palce
indices$fs=NA
indices$fs[indices$survey=="Juvenile"] = 8.5/12 
indices$fs[indices$survey=="Spring"] = 5.5/12
indices$fs[indices$survey=="Fall"] = 10.5/12

indices$log.index = NA
indices$wt = 0
ind = indices$index==0
indices$wt[!ind]=1
indices$log.index[!ind] = log(indices$index[!ind])

## create some useful objects

## all the inputs for cohort model in the list pop.dat
pop.dat = list(
  catch=catch,
  mat=mat,
  wt=wt,
  Year=sort(unique(wt.vec$Year)),
  Age=sort(unique(wt.vec$Age))
)
pop.dat$M = 0.2
pop.dat$A = length(pop.dat$Age)
pop.dat$Y = length(pop.dat$Year)
pop.dat$n = length(wt.vec$Age)
pop.dat$imap = 1:pop.dat$n  
pop.dat$age_diff = pop.dat$Age-6  
pop.dat$abs_age_diff = abs(pop.dat$Age-6)
pop.dat$veccatch = as.vector(as.matrix(pop.dat$catch)) 
pop.dat$sqrt_veccatch = sqrt(pop.dat$veccatch)

## setup the survey observation model

indices.nz = subset(indices,wt==1); #use data with estimation wt=1

temp1 = paste(wt.vec$Year,"_",wt.vec$Age,sep='')
temp2 = paste(indices.nz$Year,"_",indices.nz$Age,sep='')
indices.nz$imap = pop.dat$imap[match(temp2,temp1)]

#define q parmeters and map
indices.nz$qparm = paste(indices.nz$survey,"_",indices.nz$Age,sep='')
q.name = unique(indices.nz$qparm)
nq = length(q.name)
indices.nz$qmap = (1:nq)[match(indices.nz$qparm,q.name)]

#define index std parmseters and map
indices.nz$stdparm = indices.nz$survey; ##self-weight by survey  
#indices.nz$stdparm = 'all'
std.name = unique(indices.nz$stdparm)
nstd = length(std.name)
indices.nz$stdmap = (1:nstd)[match(indices.nz$stdparm,std.name)]

## test map
#temp=cbind(indices.nz$index,indices.nz$Year,indices.nz$Age,wt.vec$Year[indices.nz$imap],wt.vec$Age[indices.nz$imap])
#any(temp[,2]!=temp[,4]) 
#any(temp[,3]!=temp[,5])

#For now sort of cheat and take starting values from assessment

qfall = c(11,11,8,7,5,4,4,3,3)/10 
qspring = c(10,14,7,5,3,3,3,3,4)/10 
qjuvenile = c(35,18,13,11,8,6,5,3,3)/10
q.start = c(qfall,qspring,qjuvenile)      
names(q.start)=q.name

## equal weighting all indices => common log error variance parameter
std.start=c(0.3,0.3,0.3,0.1)
names(std.start) = c(std.name,'Catch')

##starting values
Nfirst = c(53067,92911,19327,16484,12049,4268,3076,3217,2287,324)/1000
No = c(63623,98989,130098,94606,135041,195489,252970,221171,121541,154111,96818,
101649,74517,42189,44126,27764,32970,54572,50070,20911,23722,33074,26374,42559,49825,
39693,10693,7819,15588,15505,6207,6865,24684,7868,801,505,969,1342,466,2809,5975,5554,
2168,990,889,1704,4661,4405,8111,13375,2704,5508,4393,1497,2174,1411,2400,5503)/1000
sparm=c(0,0) 

## spline basis functions for f_year
Xf = bs(pop.dat$Year[1:(pop.dat$Y-1)],knots=seq(1960,2015,by=3))
Xf = cbind(rep(1,nrow(Xf)),Xf)
pop.dat$Xf = Xf

fbeta = as.vector(c(log(0.2),rep(0,ncol(Xf)-1)))
f_year = exp(Xf%*%fbeta)


Npop = function(Nfirst,No,sparm,fbeta,pop.dat){ 

  A=pop.dat$A
  Y=pop.dat$Y
  M=pop.dat$M
  
#  s = exp(sparm[1]*pop.dat$age_diff + sparm[2]*pop.dat$age_diff*pop.dat$age_diff) 
  s = exp(sparm[1]*pop.dat$age_diff + sparm[2]*pop.dat$abs_age_diff)
  f_year = exp(pop.dat$Xf%*%fbeta)
  F = f_year %*% s
  Z = F + M;

  N = matrix(NA,nrow=Y,ncol=A)  
  all.age = 1:A  
  next.age = 2:A
  N[1:(Y-1),1] = No
  N[Y,1] = exp(mean(log(N[(Y-3):(Y-1),1])))
  N[1,next.age] = Nfirst
  
  for (y in seq(2,Y)){
    for (a in next.age){N[y,a] <- N[y-1,a-1]*exp(-Z[y-1,a-1])}
  }
  C = F*(1-exp(-Z))*N[1:(Y-1),]/Z 
  B = pop.dat$wt*N
  SSB = pop.dat$mat*B
  pop = list(N=N,C=C,B=B,SSB=SSB,F=F)
  return(pop)
}

create.parm = function(x){
  pname = names(x)  
  ret = list(
        Nfirst = exp(x[substr(pname,1,5)=="logNf"]),
        No=exp(x[substr(pname,1,5)=="logNo"]),      
        sparm=x[substr(pname,1,5)=="sparm"], 
        fbeta=as.vector(x[substr(pname,1,5)=="fbeta"]),
        q=exp(x[substr(pname,1,4)=="logq"]),
        std=exp(x[substr(pname,1,6)=="logstd"]))
  return(ret)
}
      
fit = function(parms,pop.dat,indices.nz){
  x = create.parm(parms) 
  pop = Npop(x$Nfirst,x$No,x$sparm,x$fbeta,pop.dat)
  N=pop$N 
  vecN = as.vector(N);
  ## make sure this works right, and years and ages are vec the same as wt.vec
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap])
  std_log_index = x$std[indices.nz$stdmap] 
  nll = -sum(dnorm(indices.nz$log.index,Elog_index,std_log_index,log=TRUE))
  nll = nll - sum(dnorm(pop.dat$sqrt_veccatch,as.vector(sqrt(pop$C)),x$std[4],log=TRUE)) 
  return(nll)
}       
      
resid = function(parms,pop.dat,indices.nz){
  x = create.parm(parms) 
  pop = Npop(x$Nfirst,x$No,x$sparm,x$fbeta,pop.dat)
  N = pop$N 
  vecN = as.vector(N);
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap])
  std_log_index = x$std[indices.nz$stdmap] 
  resid = indices.nz$log.index - Elog_index
  names(resid) = paste(indices.nz$survey,'_',indices.nz$Year,'_',indices.nz$Age,sep='')
  s.resid = resid/std_log_index
  
  Cresid = pop.dat$sqrt_veccatch-as.vector(sqrt(pop$C))  
  names(Cresid) = paste('C_',catch.vec$Year,'_',catch.vec$Age,sep='') 
  s.Cresid = Cresid/x$std[4]
  
  ret = list(resid=resid,resid_std=s.resid,Cresid=Cresid,s.Cresid=s.Cresid)
  
  #ret = data.frame(resid=c(resid,Cresid),resid_std = c(s.resid,s.Cresid))
  return(ret)
}

pred = function(parms,pop.dat,indices.nz){
  x = create.parm(parms) 
  pop = Npop(x$Nfirst,x$No,x$sparm,x$fbeta,pop.dat)
  N = pop$N  
  vecN = as.vector(N);
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap]) 
  names(Elog_index) = paste(indices.nz$survey,'_',indices.nz$Year,'_',indices.nz$Age,sep='')
  Cpred = as.vector(pop$C)  
  names(Cpred) = paste('C_',catch.vec$Year,'_',catch.vec$Age,sep='') 
  ret = list(Elog_index=Elog_index,Cpred=Cpred)
  return(ret)
}

start.parms.list = list(
  logNfirst = log(Nfirst),
  logNo = log(No),    
  sparm = sparm,
  fbeta = fbeta,
  logq = log(q.start),
  logstd = std.start
)
start.parms = unlist(start.parms.list)

lower = list(  
  logNfirst = rep(-10,length(start.parms.list$logNfirst)),
  logNo = rep(-10,length(start.parms.list$logNo)),    
  sparm = rep(-10,length(start.parms.list$sparm)),
  fbeta = rep(-10,length(start.parms.list$fbeta)),
  logq = rep(-10,length(start.parms.list$logq)),
  logstd = log(rep(0.01,length(start.parms.list$std)))
)

upper = list(  
  logNfirst = rep(100,length(start.parms.list$logNfirst)),
  logNo = rep(100,length(start.parms.list$logNo)),    
  sparm = rep(100,length(start.parms.list$sparm)),
  fbeta = rep(10,length(start.parms.list$fbeta)),
  logq = rep(10,length(start.parms.list$logq)),
  logstd = log(rep(10,length(start.parms.list$std)))
)
lower=unlist(lower)
upper=unlist(upper) 

## check initial fit
fit(start.parms,pop.dat,indices.nz)
## check gradient
temp=grad(fit,start.parms,,,,pop.dat,indices.nz)
min(abs(temp))

system.time(fit(start.parms,pop.dat,indices.nz))  

## Do the estimatiom
system.time(
model.fit <- nlminb(start.parms,fit,,,pop.dat,indices.nz,
control=list(eval.max=10000,iter.max=1000),lower=lower,upper=upper)
)

temp = grad(fit,model.fit$par,,,,pop.dat,indices.nz)   
min(abs(temp))

##get hessian
Hfit = hessian(fit,model.fit$par,,,pop.dat,indices.nz)  

parms = data.frame(log_est=model.fit$par,log_se = sqrt(diag(solve(Hfit))))
parms$est = exp(parms$log_est)
parms$se = parms$log_se*parms$est
pname = row.names(parms)  

parms$L = exp(parms$log_est - qnorm(0.975)*parms$log_se)       
parms$U = exp(parms$log_est + qnorm(0.975)*parms$log_se)

ind = (substr(pname,1,5)=='fbeta')|(substr(pname,1,5)=='sparm')  
parms$est[ind] = parms$log_est[ind]
parms$se[ind] = parms$log_se[ind] 
parms$L[ind] = parms$log_est[ind] - qnorm(0.975)*parms$log_se[ind]       
parms$U[ind] = parms$log_est[ind] + qnorm(0.975)*parms$log_se[ind]
 
pname = gsub('log','',pname)          
pname = gsub('Spring','Spr',pname)        
pname = gsub('Juvenile','Juv',pname)      
pname = gsub('Catch','C',pname) 
rownames(parms) = pname

stargazer(parms[,3:6], type = "html",  out="out.doc",align = TRUE,   
 title="3LNO cod SCA", single.row=TRUE, summary=FALSE)
# ci=TRUE, ci.level=0.95)

x = create.parm(model.fit$par)
pop = Npop(x$Nfirst,x$No,x$sparm,x$fbeta,pop.dat)


ind = match(4:6,pop.dat$Age)
aveF = apply(pop$F[,ind],1,mean)

assmt.F = read.table(file=paste(dname,'assmt_F.txt',sep=''),header=FALSE,sep="\t",col.names=c('Year','aveF'))

gname <- 'aveF.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)

  par(mar=c(3,3,0.2,0.2))
  plot(assmt.F$Year,aveF,type='l',lwd=2,xlab='',ylab='')
  lines(assmt.F$Year, unlist(assmt.F[,2]),type='l',lwd=2,col='red') 
  abline(h=0.3,col='darkorchid1',lwd=2)
  legend("topright",col=c('black','red'),lwd=2,legend=c('R SCA','NAFO ADAPT'),bty='n',cex=0.7)
  text(max(assmt.F$Year)-3,0.35,"Flim",col='darkorchid1')
  mtext(side=2,outer=F,line=2,'Ave F (Ages 4-6)')  
  mtext(side=1,outer=F,line=2,'Year') 
  
dev.off()

ssb = apply(pop$SSB,1,sum)
assmt.ssb = read.table(file=paste(dname,'assmt_ssb.txt',sep=''),header=FALSE,sep="\t",col.names=c('Year','ssb'))

ylim=range(ssb,assmt.ssb$ssb/1000)
#ylim=c(0,200)

gname <- 'ssb.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)

  par(mar=c(3,3,0.2,0.2))
  plot(assmt.ssb$Year,ssb,type='l',lwd=2,xlab='',ylab='',ylim=ylim)
  lines(assmt.ssb$Year, assmt.ssb$ssb/1000,type='l',lwd=2,col='red') 
  legend("topright",col=c('black','red'),lwd=2,legend=c('R SCA','NAFO ADAPT'),bty='n',cex=0.7)
  mtext(side=2,outer=F,line=2,'SSB (Kt)')  
  mtext(side=1,outer=F,line=2,'Year') 
  abline(h=60,col='darkorchid1',lwd=2)    
  text(max(assmt.ssb$Year)-3,65,"Blim",col='darkorchid1')
  
dev.off()

gname = "SSB_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
  persp(as.numeric(pop.dat$Year),as.numeric(pop.dat$Age),as.matrix(pop$SSB), 
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

gname = "F_3D.jpeg"
jpeg(file=gname,width=5,height=5,units='in',res=300)

  par(mar=c(0,2,1,0),mgp=c(2,0.7,0))
  persp(as.numeric(pop.dat$Year[1:(pop.dat$Y-1)]),as.numeric(pop.dat$Age),as.matrix(pop$F), 
      main="F-at-age Surface",
      xlab='Year',ylab='Age',zlab='F',
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
  persp(as.numeric(pop.dat$Year),as.numeric(pop.dat$Age),as.matrix(pop$N), 
      main="N-at-age Surface",
      xlab='Year',ylab='Age',zlab='Numbers (millions)',
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



fit.all = resid(model.fit$par,pop.dat,indices.nz)
pred.all = pred(model.fit$par,pop.dat,indices.nz)
fit = data.frame(resid = fit.all$resid,resid_std=fit.all$resid_std)

fit$Age = indices.nz$Age    
fit$Year = indices.nz$Year
fit$survey = indices.nz$survey
fit$pred = pred.all$Elog_index
fit$index = indices.nz$index

fit$colr = 'deepskyblue'
fit$colr[fit$resid>0]='firebrick2'

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year|survey, data=fit, 
       ylab='Age',main='Residuals',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 1.5*sqrt(abs(fit$resid)/pi), fill.color = fit$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(h=1:6, v=c(2000,2005,2010,2015), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
#print(resid.plot1)
                  
gname = c("resid_matrix_bubbles.jpeg")
jpeg(file=gname,width=4,height=4,units='in',res=300)
  print(resid.plot1)
dev.off()

qdat = subset(parms,substr(rownames(parms),1,1)=='q')
temp = strsplit(rownames(qdat),"_")
qdat$Age=as.numeric(unlist(lapply(temp,function(x){x[2]})))
qdat$survey=gsub('q.','',unlist(lapply(temp,function(x){x[1]})))

yttl <- list(c("Index Catchabilities (000s)"), cex=1)
xttl <- list(c("Year"), cex=1)
stripttl <- list(cex=0.2,font=1)
ax <- list(cex=0.5,relation="free")
ylim=c(0,10)

qplot = xyplot(est+L+U~Age|factor(survey), data=qdat, type="l", xlab=xttl,lwd=2,
ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(2,2,1),las=1,col=c('black','grey','grey'),
par.strip.text = list(cex=0.6),
par.settings = my.padding)

#print(qplot)
              
gname = c("q.jpeg")
jpeg(file=gname,width=3,height=3,units='in',res=300)
  print(qplot)
dev.off()

ind = wt.vec$Year<=2016
compare = data.frame(Year=wt.vec$Year[ind],Age=wt.vec$Age[ind],F=as.vector(pop$F))

  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Fishing Mortality rates (F)"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  ax <- list(cex=0.5,relation="same",tck=0.3)

jpeg(file="F_age.jpeg",width=5,height=5,units='in',res=300)
  
  print(xyplot(F~Year|factor(Age), data=compare, type="l", xlab=xttl,lwd=2,lty=1,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),las=1,
  #key=my.key,
  col=c('blue'),par.strip.text = list(cex=0.6),
  par.settings = my.padding))
  
dev.off()

age.res=compare
age.res$age=age.res$Age 
age.res$year=age.res$Year
pr.type='F';
gname1 <- "pr_F.jpeg"

bp<-function(x,y,v, scale=3, ...){
  plot(x,y,cex=sqrt(abs(v))*scale, col=ifelse(v<0,'black','grey'), pch=ifelse(v<0,16,16), ...)
  points(x[v>0],y[v>0],cex=sqrt(v[v>0])*scale, col='grey', pch=1, ...)
  points(x[v<0],y[v<0],cex=sqrt(abs(v[v<0]))*scale, col='black', pch=1, ...)
}

n.age <- length(unique(age.res$age))
n.year <- length(unique(age.res$year))

Fm <- pop$F
pr <- age.res$F/rep(tapply(age.res$F,age.res$year,max),times=n.age)

Prm <- matrix(pr,ncol=n.age,byrow=FALSE) 
ave.pr <- tapply(pr,age.res$age,mean)

ave.Prm <- matrix(ave.pr,nrow=n.year,ncol=n.age,byrow=TRUE)
Prdiff <- Prm - ave.Prm;
                                
yval <- unique(age.res$age) 
xval <- unique(age.res$year) 

jpeg(file=gname1,width=5,height=6,units='in',res=300)

par(mfrow=c(2,1),oma=c(1,2,0,0),mar=c(3,1,1,0))
nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(7), c(3,6), TRUE)

plot(yval,ave.pr/max(ave.pr),xlab='',ylab='',type='l',lwd=2,las=1,cex=0.8)
mtext(side=3,outer=F,line=0.1,'Average Selectivity Pattern',cex=0.8)
mtext(side=1,outer=F,line=1.8,'Age',cex=1)  

vPrdiff <- as.vector(Prdiff); # this is in age within year order;

scale <- 2
 bp(age.res$year,age.res$age,vPrdiff, ylim=range(age.res$age), xlim=range(age.res$year), 
       las=1, xlab='', ylab='', main="", scale=scale)
    points(age.res$year,age.res$age,cex=sqrt(abs(vPrdiff))*scale, col='black', pch=1)
mtext(side=2,outer=F,line=2,'Age',cex=1) 
mtext(side=1,outer=F,line=2.5,'Year',cex=1) 
mtext(side=3,outer=F,line=0.1,'Change in Selectivity Pattern',cex=0.8)

dev.off()

sname = unique(fit$survey)
fit.all=fit
plot_type='jpeg'

fit = fit.all
fit$year=fit$Year
fit$age=fit$Age
fit$Eindex = exp(fit$pred)

tran.bias = 1
gname <- 'index.jpeg'
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\obs_pred_plot.txt",local=FALSE)
dev.off()


fit=fit.all 
fit$year=fit$Year
fit$age=fit$Age
fit$Eindex = exp(fit$pred)
usurvey=sname
fit$std.resid = fit$resid_std
fit$cohort=fit$year-fit$age
fit$est.wt=1
uage = sort(unique(fit$age))
pname='';fig.folder=getwd()
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\sres_plot.txt",local=FALSE)

fit$log_index = indices.nz$log.index
fit$Elog_index = log(fit$Eindex)
   
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\fit_plot.txt",local=FALSE)

ht <- 6; # height of plot in inches;
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\sres_smooth_plot.txt",local=FALSE)

#################  catch fit plots #####################;

fit.all = resid(model.fit$par,pop.dat,indices.nz)
pred.all = pred(model.fit$par,pop.dat,indices.nz)

fit = data.frame(resid = fit.all$Cresid,resid_std=fit.all$s.Cresid)

fit$Age = catch.vec$Age    
fit$Year = catch.vec$Year
fit$pred = pred.all$Cpred
fit$index = catch.vec$index

fit$colr = 'deepskyblue'
fit$colr[fit$resid>0]='firebrick2'

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year, data=fit, 
       ylab='Age',main='Catch Standardized Residuals',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 1.5*sqrt(abs(fit$resid)/pi), fill.color = fit$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(h=1:6, v=seq(1960,2015,by=5), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
#print(resid.plot1)
                  
gname = c("catch_resid_matrix_bubbles.jpeg")
jpeg(file=gname,width=4,height=4,units='in',res=300)
  print(resid.plot1)
dev.off()



catch.tot =  aggregate(cbind(fit$index,fit$pred),list(Year=fit$Year),sum)
colnames(catch.tot) = c('Year','O','P')

ylim=range(catch.tot[,2:3])

gname <- 'obs_pred_total_catch.jpeg'
jpeg(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,0.5,0.5),cex.axis=1,las=1,mgp=c(1,0.5,0))

plot(catch.tot$Year,catch.tot$P,type='l',las=1,xlab='',ylab='',lwd=2,ylim=ylim)
points(catch.tot$Year,catch.tot$O,pch=16)       
lines(catch.tot$Year,catch.tot$P,lty=1,lwd=2,col='red')
 
yname = c("Total Catch (millions)")
mtext(side=2,line=2,yname,las=0,cex=1)  
mtext(side=1,line=2,'Year',las=0,cex=1)


legend('topright',bty='n',c('Observed','Predicted'),lty=c(NA,1),lwd=c(NA,2),pch=c(16,NA),
col=c('black','red'))

dev.off()

ttl <- list(c("Catch at ages 2-12"), cex=1)
  xttl <- list(c("Year"), cex=1)
  yttl <- list(c("Catch (millions)"), cex=1)
  stripttl <- list(cex=0.2,font=1)
  ax <- list(cex=0.5,relation="free")
  my.key <- list(space = "top",columns=2,between=1,
                border = FALSE,
                size = 3,
                lines = list(lty = 1,col=2:1,lwd=2),
                text = list(c("Observed","Predicted")))

gname <- 'obs_pred_catch.jpeg'
jpeg(file=gname,width=4,height=5,units='in',res=300)

  print(xyplot(index+pred~Year|factor(Age), data=fit, type="l", xlab=xttl,lwd=2,
  ylab=yttl, strip.text=stripttl, scales=ax, as.table=TRUE, layout=c(3,4,1),las=1,
  key=my.key,col=2:1,par.strip.text = list(cex=0.6),
  par.settings = list(layout.heights = list(strip = .8))))

dev.off()

fit$year=fit$Year
fit$age=fit$Age
fit$survey='Catch'
fit$Eindex = fit$pred
usurvey='Catch'
fit$std.resid = fit$resid_std
fit$cohort=fit$year-fit$age
fit$est.wt=1
ind = fit$index==0
fit$est.wt[ind]=0
uage = sort(unique(fit$age))
pname='';fig.folder=getwd()
source("C:\\home\\CADIGAN\\Rfunc\\admb_plots\\sres_plot.txt",local=FALSE)

