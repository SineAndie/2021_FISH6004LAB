# lab 6

# lognormal VonB data
n = 1000
age = sort(rpois(n,lambda=6))
barplot(table(age),xlab='Age',ylab='Frequency',main='Age Distribution')
Linf = 55
K = 0.35
ao = -.1
sigma=0.2
mu = Linf*(1 - exp(-K*(age - ao)))
len = rlnorm(n, meanlog = log(mu), sdlog = sigma)
plot(age,len,xlab='Age',ylab="Length",type='p',xlim=c(0,20))
lines(age,mu,lty=1,lwd=2,col='red')

vonb <- function(Linf,K,ao,age){
  mu = Linf*(1 - exp(-K*(age-ao)))
  return(mu)
}

lognormal.nll <- function(parm){
  Linf = exp(parm[1])
  K = exp(parm[2])
  ao = parm[3]
  sigma = exp(parm[4])
  mu = vonb(Linf,K,ao,age)
  prob = dlnorm(len,mean=log(mu),sd=sigma)
  nll = -sum(log(prob))
  return(nll)
}

parm.start = c(3,-1,-0.1,-1)
mle.fit = nlminb(parm.start,lognormal.nll)

#plot model prediction 
pred.age = seq(0,20,by=0.1)

Linf.hat = exp(mle.fit$par[1])
K.hat = exp(mle.fit$par[2])
ao.hat = mle.fit$par[3]

pred.mu = vonb(Linf.hat,K.hat,ao.hat,pred.age)

lines(pred.age,pred.mu,lty=1,lwd=2,col='blue')
legend("bottomright",c("True","Estimated"),
       lty=1,col=c('red','blue'),lwd=2)

# result 
print(age)
print(len)
parm.start = c(3,-1,-0.1,-1)
mle.fit = nlminb(parm.start,lognormal.nll)
print(mle.fit$par)
print(mle.fit$converge)
vonb.parms = c(exp(mle.fit$par[1]),exp(mle.fit$par[2]),
               + mle.fit$par[3])
names(vonb.parms) = c('Linf','K','ao')
print(vonb.parms,digits=2)
sigma = exp(mle.fit$par[4])
print(sigma)

# Gamma vonB
gamma.nll <- function(parm){
  Linf = exp(parm[1])
  K = exp(parm[2])
  ao = parm[3]
  sigma = exp(parm[4])
  mu = vonb(Linf,K,ao,age)
  var = (sigma*mu)**2
  s = var/mu
  a = mu/s;
  prob = dgamma(len,shape=a,scale=s)
  nll = -sum(log(prob))
  return(nll)
}

parm.start = c(5,-1,-0.1,-1)
mle.fit = nlminb(parm.start,gamma.nll)
print(mle.fit$par)
print(mle.fit$converge)
vonb.parms = c(exp(mle.fit$par[1]),exp(mle.fit$par[2]),
               + mle.fit$par[3])
names(vonb.parms) = c('Linf','K','ao')
print(vonb.parms,digits=2)
sigma = exp(mle.fit$par[4])
print(sigma)

# Binomial
matdata <- read.csv("fname.csv",header = TRUE)
print(matdata)
plot(matdata$age,matdata$pro,xlab='Age',ylab="Proportion Mature",type='p')

pmat <- function(parm,age){
  logitp = parm[1] + parm[2]*age
  prob.mat = exp(logitp)/(1+exp(logitp))
  return(prob.mat)
}
binomial.nll <- function(parm){
  prob.mat = pmat(parm,age);
  prob = dbinom(y,n,prob.mat)
  nll = -sum(log(prob))
  return(nll)
}
y = matdata$y
n = matdata$n
age = matdata$age

parm.start = c(-4,1)
mle.fit = nlminb(parm.start,binomial.nll)
print(mle.fit$par)
print(mle.fit$converge)

pred.prob =
  pmat(mle.fit$par,age)
lines(matdata$age,pred.prob
      ,type='l',lwd=2,col='red')
age.pred =
  seq(min(matdata$age),max(matdata$age),length=200)
pred.prob1 = pmat(mle.fit$par,age.pred)
lines(age.pred,pred.prob1,type='l',lwd=2,col='blue')
legend('bottomright',lty=1,lwd=2,col=c('black','red','blue'),
       c("observed","predicted","smooth predicted"))

# Binomial logistic regression
Binfit <- glm(cbind(y,n-y) ~
                age,family=binomial(link = "logit"),data=matdata)
temp = summary(Binfit)
print(temp$coefficients)
confint.default(Binfit)
confint(Binfit)

# bootstrap
mu=10;
sigma2=1
n=10
y = rnorm(n,mu,sd=sqrt(sigma2))
ybar=mean(y)
s2 = var(y)

B=5000
booty = lapply(1:B, function(i) sample(y, replace = T))
boot.mean = sapply(booty,mean)
boot.var = sapply(booty,var)
boot.t = (boot.mean - ybar)/sqrt(boot.var/n)
quant.boot.t = quantile(boot.t,probs=c(0.025,0.975))

hist(boot.t,breaks=50,xlim=c(-4,4),freq=F)
x=seq(-4,4,length=100)
pt = dt(x,n-1)
lines(x,pt,lwd=2,col='red',lty=1)

hist(boot.mean,breaks=50,freq=F)
x=seq(8,11,length=100)
pz = dnorm(x,mean=ybar,sd=sqrt(s2/n))
lines(x,pz,lwd=2,col='red',lty=1)





