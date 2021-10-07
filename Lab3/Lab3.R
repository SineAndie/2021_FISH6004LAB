#set working directory
setwd("")

#need this package to calculate derivs
require("numDeriv")

########################################################PART I
#1) Simulate logistic population dynamics
# initial population size and parameters r, K
n0 = 1
r = 0.5
K = 100

logistic.continous <- function (n0, r, K, T) {

  
}

T = seq(0,100, length.out = 1000)
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(0, 100), 
     ylim = c(0, 100))

n1 = logistic.continous(n0,r,K,T)

lines(y=n1,x=T,col="blue")

n.change = diff(n1)
plot(n1[1: length(n1)-1], n.change)
abline(v = K/2, col = "green")

#try for r = 1.5, 2.9


#2) Surplus production model
b0=1
r=1
K=100
T=100
H=0

surplus.production <- function (b0, r, K, T, H) {

  
}

b = surplus.production(b0,r,K,T,H)
par(mfrow=c(1,2))
plot(b)
plot(b[1: length(b)-1], b[2: length(b)], xlab = "Bt", ylab ="B(t+1)")
abline(0, 1, col = "red") 
abline(v = K/2, col = "green")

# grid out B, find Fmsy
B = seq(0, 100, length.out = 1000)
C = r*B*(1-B/K)
plot(B,C)

F = r*(1-B/K)
plot(F,C)
F[which.max(C)]



#3) 
############### YPR ##################################;
# uniform selectivity
YPR <- function(f){
  n.age = length(sel)
  fa = f*sel
  za = m + fa
  cza = cumsum(c(0,za[1:(n.age-1)]))
  YPR = sum(exp(-cza)*weight*fa*(1-exp(-za))/za)
  return(YPR)
}

age=3:20
sel = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

weight = c(0.15,0.278,0.394,0.599,0.862,1.163,1.572,2.028,2.653,
           3.141,3.844,4.043,4.252,4.471,4.702,4.945,5.200,5.469)

m=rep(0.2,length(age))

ypr1 = YPR(0.2)
print(ypr1)

## YPR for a range of F's; change sel and see what happens to Fmax
fv = seq(0,0.5,by=0.01)
fv = as.matrix(fv,ncol=1)
ypr2 = apply(fv,1,YPR)

plot(fv,ypr2,xlab='F',ylab='YPR',type='l',ylim=c(0,0.4))

Fmax = fv[which.max(ypr2)]
abline(v=Fmax)
text(Fmax+0.03,0.1,'Fmax')

sel = 
ypr3 = apply(fv,1,YPR)
lines(fv,ypr3,lty=1,col='red')
Fmax = fv[which.max(ypr3)]
abline(v=Fmax,col='red')
text(Fmax+0.03,0.12,'Fmax',col='red')

sel = 
ypr3 = apply(fv,1,YPR)
lines(fv,ypr3,lty=1,col='green')
Fmax = fv[which.max(ypr3)]
abline(v=Fmax,col='green')
text(Fmax+0.03,0.14,'Fmax',col='green')

## try changing other inputs (e.g. m, weights) and assess what happens to Fmax

#now calculate F0.1
#####   F0.1  ##################
fv = seq(0,0.5,by=0.001)
fv = as.matrix(fv,ncol=1)
ypr = apply(fv,1,YPR)

plot(fv,ypr,xlab='F',ylab='YPR',type='l')

Fmax = fv[which.max(ypr)]
abline(v=Fmax)
text(Fmax+0.03,0.1,'Fmax')

der.ypr = rep(0,length(ypr))
for(i in 1:length(ypr)){
  der.ypr[i] = grad(YPR,fv[i])
}
der.ypr.scaled = der.ypr/der.ypr[1]
ipos = which.min(abs(der.ypr.scaled-0.1))
F0.1 = fv[ipos]
abline(v=F0.1,lty=2)
text(F0.1+0.03,0.1,'F0.1')

YPR0.1 = ypr[ipos]
b = der.ypr[ipos]; ## slope at F0.1
a = YPR0.1 - b*F0.1  ## intercept of tangent line
abline(a,b,lty=1,col=grey(0.7))

YPR1 = ypr[1]
b = der.ypr[1]; ## slope at origin
a = YPR1  ## intercept of tangent line
abline(a,b,lty=2,col=grey(0.7))

ctext = paste('slope at origin = ',signif(der.ypr[1],digits=2),
              ', slope at F0.1 = ',signif(der.ypr[ipos],digits=2),sep='')
mtext(side=3,line=0.5,ctext)

###################### SPR #############################;
SPR <- function(f){
  n.age = length(sel)
  fa = f*sel
  za = m + fa
  cza = cumsum(c(0,za[1:(n.age-1)]))
  SPR = sum(exp(-cza)*weight*mat)
  return(SPR)
}

age=3:20
sel = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

weight = c(0.15,0.278,0.394,0.599,0.862,1.163,1.572,2.028,2.653,
           3.141,3.844,4.043,4.252,4.471,4.702,4.945,5.200,5.469)

mat = c(0,0,0,0,0,0.002,0.006,0.022,0.077,0.242,0.529,0.651,0.756,
        0.837,0.895,0.934,0.959,0.975)

m=rep(0.2,length(age))

spr1 = SPR(0.2)
print(spr1)

## SPR for a range of F's; change sel and see what happens 
fv = seq(0,0.5,by=0.001)
fv = as.matrix(fv,ncol=1)
spr1 = apply(fv,1,SPR)
percent.spr1 = 100*spr1/spr1[1]

plot(fv,percent.spr1,xlab='F',ylab='%SPR',type='l',ylim=c(0,100))

ipos=which.min(abs(percent.spr1-35))
F35 = fv[ipos]
segments(F35,0,F35,35)
text(F35+0.03,35,'F35',cex=0.8)
ipos=which.min(abs(percent.spr1-20))
F20 = fv[ipos]
segments(F20,0,F20,20)
text(F20+0.03,20,'F20',cex=0.8)

sel = 
spr2 = apply(fv,1,SPR)
percent.spr2 = 100*spr2/spr2[1]
lines(fv,percent.spr2,lty=1,col='green')

ipos=which.min(abs(percent.spr2-35))
F35 = fv[ipos]
segments(F35,0,F35,35,col='green')
text(F35+0.03,35,'F35',cex=0.8,col='green')
ipos=which.min(abs(percent.spr2-20))
F20 = fv[ipos]
segments(F20,0,F20,20,col='green')
text(F20+0.03,20,'F20',cex=0.8,col='green')

ipos=which.min(abs(fv-0.3))
abline(h=percent.spr2[ipos],lty=2,col='green')

#4##############  Beverton-Holt and MSY ######################;

BH = function(alpha,beta,ssb){
  alpha*ssb/(beta+ssb)
}

alpha=1000
beta=10
ssb=seq(1:200)

rec = BH(alpha,beta,ssb)

plot(ssb,rec,type='l',lwd=2)

## can show for BH that equilibrium SSB for a constant fishing mortality f is (see slides);

SSBe = function(f){alpha*SPR(f)-beta}

RECe = function(f){BH(alpha,beta,SSBe(f))}

Yielde = function(f){RECe(f)*YPR(f)}

fv = seq(0,0.5,by=0.001)
fv = as.matrix(fv,ncol=1)
ye_f = apply(fv,1,Yielde)

plot(fv,ye_f,xlab='F',ylab='Equilibrium Yield',type='l')

#Calculate Fmsy
Fmsy = 

#calculate Bmsy
Bmsy = 

#------------------------------------------------------------------
## Continuous time YPR;

## The Von Bertalanffy weight at age t function;
vonb = function(t){
  ret = Winf*((1 - exp(-k*t))**b)
  return(ret)
}  

ypr.int = function(t,f){
  
  wt = vonb(t)  
  st = pnorm(t-tc); ## selectivity as a normal pdf;
  St = (t-tc)*pnorm(t-tc) + dnorm(t-tc) + tc*pnorm(-tc) - dnorm(-tc)
  Mt = m*t; ## m(t)=m so that M(t) = int_0^t m(t) dt = m*t                       
  ft = f*St
  ret = f*st*wt*exp(-ft-Mt);
  return(ret)                                      
}

ypr.fun = function(f){integrate(ypr.int,f, lower = 0, upper = Inf)[[1]]}

b=3
Winf = 25
k = 0.2 
m=0.2
tc=4

fv = seq(0.1,1,by=0.01)
fv = as.matrix(fv,ncol=1)
YPRf = apply(fv,1,ypr.fun)

plot(fv,YPRf,type='l') 

Fmax = fv[which.max(YPRf)]
abline(v=Fmax)
text(Fmax+0.03,0.1,'Fmax')

## Another way to get Fmax, define a function that returns dYPR(f)/df

dyprt.df = function(t,f){  
  st = pnorm(t-tc)
  St = (t-tc)*pnorm(t-tc) + dnorm(t-tc) + tc*pnorm(-tc) - dnorm(-tc)          
  ft = f*St
  ret = ypr.int(t,f)*(1-ft)/f;
  return(ret)
}
upper=Inf
dypr.df = function(f){integrate(dyprt.df,f, lower = 0, upper = upper)[[1]]}  

Fmax1 = uniroot(dypr.df, interval=c(1e-8,1),tol=1e-8)$root

abline(v=Fmax1,col='red')

### SPR  ####

spr.int = function(t,f){
  wt = vonb(t)
  st = pnorm(t-tc)
  St = (t-tc)*pnorm(t-tc) + dnorm(t-tc) + tc*pnorm(-tc) - dnorm(-tc)
  Mt = m*t                       
  ft = f*St
  logit_pt =  bo+b1*t
  pt=rep(1,length(t)); ## proportion mature a logistic function of t;
  ind = logit_pt<100
  pt1= exp(bo+b1*t[ind])/(1+exp(bo+b1*t[ind]))
  pt = replace(pt,ind,pt1)
  ret = pt*wt*exp(-ft-Mt);
  return(ret)
} 

spr.fun = function(f){integrate(spr.int, lower = 0, upper = Inf,f)[[1]]}

dsprt.df = function(t,f){
  st = pnorm(t-tc)
  St = (t-tc)*pnorm(t-tc) + dnorm(t-tc) + tc*pnorm(-tc) - dnorm(-tc) 
  ret = -spr.int(t,f)*St;
  return(ret)
} 


t50=4; #Age at 50% maturity;
t95=5.5; #age at 95% maturity;
b1=log(0.95/0.05)/(t95-t50)
bo = -t50*b1 


## compute the percent SPR, 100*SPR(f)/SPR(0)
SPRf = 100*apply(fv,1,spr.fun)/spr.fun(0)

plot(fv,SPRf,type='l') 

ipos=which.min(abs(SPRf-35))
F35 = fv[ipos]
segments(F35,0,F35,35)
text(F35+0.03,35,'F35',cex=0.8)
ipos=which.min(abs(SPRf-20))
F20 = fv[ipos]
segments(F20,0,F20,20)
text(F20+0.03,20,'F20',cex=0.8)
abline(h=c(20,35),lty=2)


