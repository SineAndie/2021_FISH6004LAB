
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

load('data\\3LNO.RData')

#> ls()
# "catch"     "catch.vec" "indices"   "mat"       "mat.vec"   "wt"        "wt.vec" 

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
##indices.nz$stdparm = indices.nz$survey; ##self-weight by survey  
indices.nz$stdparm = 'all'
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
std.start=0.5

##set survivors
N_Y = c(5.4,1.7,0.7,0.7,0.3,1,1,0.4,1.9,1)
names(N_Y) = paste('N_2017,',3:12,sep='')
N_Aest = c(12,6,4,8,11,5,6,52,24,3,4,7,8,2,12,27,19,5,3,6,13,41,45)/100
Years_NAest = 1994:2016      
names(N_Aest) = paste('N_',Years_NAest,',12',sep='')
pop.dat$NA.map = match(Years_NAest,pop.dat$Year)
pop.dat$Fave.map = match(6:9,pop.dat$Age)

Npop = function(N_Y,N_Aest,pop.dat){ 

  A=pop.dat$A
  Y=pop.dat$Y
  C=pop.dat$catch 
  M=pop.dat$M

  N = matrix(NA,nrow=Y,ncol=A)      
  F = matrix(NA,nrow=Y-1,ncol=A) 
  last.age = 1:(A-1)           
  all.age = 1:A
  N[Y,2:A]=N_Y
  N[pop.dat$NA.map,A] = N_Aest
  
  for (y in seq(Y-1,1,-1)){
    for (i in seq_along(last.age)){
      N[y,i] <- N[y+1,i+1]*exp(M) + C[y,i]*exp(M/2)
      F[y,i] = -log(1 - exp(M/2)*C[y,i]/N[y,i]);
    }
    test = is.na(N[y,A])
    F[y,A] = ifelse(test,mean(F[y,pop.dat$Fave.map]),-log(1 - exp(M/2)*C[y,A]/N[y,A]))
    if(test){N[y,A] = C[y,A]*exp(M/2)/(1-exp(-F[y,A]))}
  } 
  N[Y,1] = exp(mean(log(N[(Y-3):(Y-1),1])))
  B = pop.dat$wt*N
  SSB = pop.dat$mat*B
  pop = list(N=N,B=B,SSB=SSB,F=F)
  return(pop)
}

create.parm = function(x){
  pname = names(x)  
  ret = list(
        N_Y = exp(x[substr(pname,1,6)=="logN_Y"]),
        N_Aest=exp(x[substr(pname,1,6)=="logN_A"]),
        q=exp(x[substr(pname,1,4)=="logq"]),
        std_index=exp(x[substr(pname,1,9)=="std_index"]))
  return(ret)
}
      
fit = function(parms,pop.dat,indices.nz){
  x = create.parm(parms)
  N = Npop(x$N_Y,x$N_Aest,pop.dat)$N 
  vecN = as.vector(N);
  ## make sure this works right, and years and ages are vec the same as wt.vec
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap])
  std_log_index = x$std_index[indices.nz$stdmap] 
  nll = -sum(dnorm(indices.nz$log.index,Elog_index,std_log_index,log=TRUE))
  return(nll)
}       
      
resid = function(parms,pop.dat,indices.nz){
  x = create.parm(parms)
  N = Npop(x$N_Y,x$N_Aest,pop.dat)$N 
  vecN = as.vector(N);
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap])
  std_log_index = x$std_index[indices.nz$stdmap] 
  resid = indices.nz$log.index - Elog_index
  s.resid = resid/std_log_index
  ret = data.frame(resid=resid,resid_std = s.resid)
  return(ret)
}

pred = function(parms,pop.dat,indices.nz){
  x = create.parm(parms)
  N = Npop(x$N_Y,x$N_Aest,pop.dat)$N 
  vecN = as.vector(N);
  Elog_index = log(x$q[indices.nz$qmap]) + log(vecN[indices.nz$imap])
  return(Elog_index)
}

start.parms.list = list(
  logN_Y = log(N_Y),
  logN_Aest = log(N_Aest),
  logq = log(q.start),
  std_index = log(0.5)
)
start.parms = unlist(start.parms.list)

lower = list(    
  logN_Y = rep(-10,length(start.parms.list$logN_Y)),
  logN_Aest = rep(-10,length(start.parms.list$logN_Aest)),
  logq = rep(-10,length(start.parms.list$logq)),
  std_index = log(0.1)
)

upper = list(   
  logN_Y = rep(100,length(start.parms.list$logN_Y)),
  logN_Aest = rep(100,length(start.parms.list$logN_Aest)),
  logq = rep(100,length(start.parms.list$logq)),
  std_index = log(10)
)
lower=unlist(lower)
upper=unlist(upper) 

## check initial fit
fit(start.parms,pop.dat,indices.nz)
## check gradient
grad(fit,start.parms,,,,pop.dat,indices.nz)

system.time(fit(start.parms,pop.dat,indices.nz))  

## Do the estimation
adapt.fit <- nlminb(start.parms,fit,,,pop.dat,indices.nz,lower=lower,upper=upper)
##get hessian
Hfit = hessian(fit,adapt.fit$par,,,pop.dat,indices.nz)  

parms = data.frame(log_est=adapt.fit$par,log_se = sqrt(diag(solve(Hfit))))
parms$est = exp(parms$log_est)
parms$se = parms$log_se*parms$est
parms$L = exp(parms$log_est - qnorm(0.975)*parms$log_se)       
parms$U = exp(parms$log_est + qnorm(0.975)*parms$log_se)
pname = gsub('log','',rownames(parms))    
pname = gsub('.N','',pname) 
ind = substr(pname,1,6)=="N_Aest" 
pname[ind] = gsub(',12','',pname[ind])       
pname = gsub('_2017','',pname)          
pname = gsub('Aest','A',pname)          
pname = gsub('Spring','Spr',pname)        
pname = gsub('Juvenile','Juv',pname) 
rownames(parms) = pname

stargazer(parms[,3:6], type = "html",  out="out.doc",align = TRUE,   
 title="3LNO cod ADAPT", single.row=TRUE, summary=FALSE)
# ci=TRUE, ci.level=0.95)

x = create.parm(adapt.fit$par)
pop = Npop(x$N_Y,x$N_Aest,pop.dat)


ind = match(4:6,pop.dat$Age)
aveF = apply(pop$F[,ind],1,mean)

assmt.F = read.table(file='data\\assmt_F.txt',header=FALSE,sep="\t",col.names=c('Year','aveF'))

gname <- 'aveF.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)

  par(mar=c(3,3,0.2,0.2))
  plot(assmt.F$Year,aveF,type='l',lwd=2,xlab='',ylab='')
  lines(assmt.F$Year, unlist(assmt.F[,2]),type='l',lwd=2,col='red') 
  abline(h=0.3,col='darkorchid1',lwd=2)
  legend("topright",col=c('black','red'),lwd=2,legend=c('R ADAPT','NAFO ADAPT'),bty='n',cex=0.7)
  text(max(assmt.F$Year)-3,0.35,"Flim",col='darkorchid1')
  mtext(side=2,outer=F,line=2,'Ave F (Ages 4-6)')  
  mtext(side=1,outer=F,line=2,'Year') 
  
dev.off()

ssb = apply(pop$SSB,1,sum)
assmt.ssb = read.table(file='data\\assmt_ssb.txt',header=FALSE,sep="\t",col.names=c('Year','ssb'))

ylim=range(ssb,assmt.ssb$ssb/1000)

gname <- 'ssb.jpeg'
jpeg(file=gname,width=4,height=3,units='in',res=300)

  par(mar=c(3,3,0.2,0.2))
  plot(assmt.ssb$Year,ssb,type='l',lwd=2,xlab='',ylab='')
  lines(assmt.ssb$Year, assmt.ssb$ssb/1000,type='l',lwd=2,col='red') 
  legend("topright",col=c('black','red'),lwd=2,legend=c('R ADAPT','NAFO ADAPT'),bty='n',cex=0.7)
  mtext(side=2,outer=F,line=2,'SSB (Kt)')  
  mtext(side=1,outer=F,line=2,'Year') 
  abline(h=60,col='darkorchid1',lwd=2)    
  text(max(assmt.ssb$Year)-3,65,"Blim",col='darkorchid1')
  
dev.off()


ssb.matrix = matrix(pred.dat$vSSB,nrow=Y,ncol=A,byrow=FALSE)

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



fit = resid(adapt.fit$par,pop.dat,indices.nz)
fit$Age = indices.nz$Age    
fit$Year = indices.nz$Year
fit$survey = indices.nz$survey
fit$pred = pred(adapt.fit$par,pop.dat,indices.nz)
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
xttl <- list(c("Age"), cex=1)
stripttl <- list(cex=0.2,font=1)
ax <- list(cex=0.5,relation="free")

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
#  my.key <- list(space = "top",columns=3,
#                border = FALSE,
#                size = 3,between=0.1,
#                lines = list(lty = c(1,1,1),col=c('red','green','blue'),lwd=c(2,2,2)),
#                text = list(c("F","Mo","M")))
  
#  ax <- list(cex=0.5,relation="free",tck=0.3)
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

Prm <- matrix(pr,ncol=n.age,byrow=F) 
ave.pr <- tapply(pr,age.res$age,mean)

ave.Prm <- matrix(ave.pr,nrow=n.year,ncol=n.age,byrow=T)
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
source("obs_pred_plot.txt",local=F)
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
source("sres_plot.txt",local=F)

fit$log_index = indices.nz$log.index
fit$Elog_index = log(fit$Eindex)
   
source("fit_plot.txt",local=F)

ht <- 6; # height of plot in inches;
source("sres_smooth_plot.txt",local=F)




