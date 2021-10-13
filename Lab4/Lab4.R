#set working directory
setwd("")

#need this package to calculate derivs
require("numDeriv")

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

##############change selectivity 
sel = c(0.2,0.5,0.8,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
ypr3 = apply(fv,1,YPR)
lines(fv,ypr3,lty=1,col='red')
Fmax = fv[which.max(ypr3)]
abline(v=Fmax,col='red')
text(Fmax+0.03,0.12,'Fmax',col='red')

sel = c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
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

## SPR for a range of F's; change sel and see what happens to Fmax
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

sel = c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
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
Fmsy = fv[which.max(ye_f)]
abline(v=Fmsy)
text(Fmsy+0.03,0.1,'Fmsy')

#calculate Bmsy
Bmsy = SSBe(Fmsy)




