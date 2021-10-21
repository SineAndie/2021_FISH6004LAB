#------------------------Part I----------------------------
# try nlmbin
#generate a sample from known population parameters
n=20
pop.mu=5
pop.s2=1;

y = rnorm(n=n,mean=pop.mu,sd=sqrt(pop.s2))
print(y,digits=3)

ybar=mean(y)
print(ybar)

yvar = sum((y-ybar)^2)/n
print(yvar)

#pretend we don't know the pop mu and s2 and estimate 
normal.nll <- function(parm){
  mu  = parm[1]
  s2  = exp(parm[2])
  nll =  n*log(2*pi*s2)/2 + sum((y-mean(y))^2)/(2*s2) + n*(mean(y)-mu)^2/(2*s2)
  return(nll)
}

parm.start = c(0,1)
normal.nll(parm.start)

mle.fit = nlminb(parm.start,normal.nll)
print(mle.fit)

exp(mle.fit$par[2])

###Same as:
normal.nll2 <- function(parm){
  mu  = parm[1]
  s2  = exp(parm[2])
  prob= dnorm(y,mean=mu,sd=sqrt(s2))
  nll = -sum(log(prob))
  return(nll)
}

parm.start = c(0,1)
normal.nll2(parm.start)
mle.fit2 = nlminb(parm.start,normal.nll2)
print(mle.fit2)

###Same as:
normal.nll3 <- function(parm){
  mu  = parm[1]
  s2  = exp(parm[2])
  nll= -sum(dnorm(y,mean=mu,sd=sqrt(s2),log = TRUE))
  return(nll)
}

parm.start = c(0,1)
normal.nll3(parm.start)
mle.fit3 = nlminb(parm.start,normal.nll3)

mle.fit$par
mle.fit2$par
mle.fit3$par

#------------------------Part II----------------------------
# code up a optimizer to understand how nlminb works
set.seed(1234)

library(numDeriv)

vonb = function(Linf,K,ao,age){
  len = Linf*(1-exp(-K*(age-ao)));
  return(len)
}

tdata=data.frame(age=seq(1,20,by=1))
n = length(tdata$age)

pop.Linf = 120
pop.K = 0.15
pop.a0 = 0.5
err.mu = 0
err.s2 = 5

tdata$len = vonb(pop.Linf,pop.K,pop.a0,tdata$age) + rnorm(n,err.mu,sd=err.s2)

plot(tdata$age,tdata$len)

vonb_nll = function(parm){
  
  Linf = exp(parm[1])
  K = exp(parm[2])
  ao = exp(parm[3])
  s2 = exp(parm[4])

  mu = vonb(Linf,K,ao,tdata$age)
  n = length(mu)
  resid = tdata$len-mu
  nll = n*log(2*pi*s2)/2 + sum(resid**2)/(2*s2)
  return(nll)
}

parm.true = log(c(pop.Linf,pop.K,pop.a0,err.s2))
nlminb(parm.true, vonb_nll)

## jacobian: matrix of all first order partial derivatives 
der.vonb = function(Linf,K,ao,age){
  mu = Linf*(1-exp(-K*(age-ao))) # mu is always the E(L)
  dmu.dLinf = 1-exp(-K*(age-ao))
  dmu.dK = Linf*exp(-K*(age-ao))*(age-ao)
  dmu.dao = -Linf*exp(-K*(age-ao))*K
  
  ret = cbind(dmu.dLinf,dmu.dK,dmu.dao)
  return(ret)
}

##check
der.vonb(120,0.1,0.1,tdata$age)


vonb_wrap =  function(parm,age){
  return(vonb(parm[1],parm[2],parm[3],age))
}

jacobian(vonb_wrap, c(120,0.1,0.1),age=tdata$age) 

## gradient: vector of first order derivatives
der.vonb_nll = function(parm){
  Linf = exp(parm[1])
  K = exp(parm[2])
  ao = exp(parm[3])
  s2 = exp(parm[4])
  mu = vonb(Linf,K,ao,tdata$age)
  n = length(mu)
  resid = tdata$len-mu
  nll = n*log(2*pi*s2)/2 + sum(resid**2)/(2*s2) 
  jac = der.vonb(Linf,K,ao,tdata$age)   
  dnll.dreg_parm = apply(-resid*jac/s2,2,sum)
  dnll.ds2 = n/(2*s2) - 0.5*sum(resid**2)/(s2**2) 
  
  ret = c(dnll.dreg_parm,dnll.ds2)
  names(ret) = c('dnll.dLinf','dnll.dK','dnll.dao','dnll.ds2') 
  
  return(ret*exp(parm))
}
# note: last step because ret is deriv with respect to natural scale parms;
# but we need derivs wrt to log(parm)

## check
start.parm = log(c(120,0.15,0.5,25))
der.vonb_nll(start.parm)

grad(vonb_nll,start.parm)

## hessian
## next three functions are derivatives of jacobian wrt to Linf, k, and ao
der.vonb.dLinf = function(Linf,K,ao,age){
  ret = cbind(0,exp(-K*(age-ao))*(age-ao),-exp(-K*(age-ao))*K)
  return(ret)
}

der.vonb.dK = function(Linf,K,ao,age){
  ret = cbind(exp(-K*(age-ao))*(age-ao),-Linf*exp(-K*(age-ao))*((age-ao)**2),Linf*exp(-K*(age-ao))*K*(age-ao)-Linf*exp(-K*(age-ao)))
  return(ret)
}

der.vonb.dao = function(Linf,K,ao,age){
  ret = cbind(-exp(-K*(age-ao))*K,Linf*exp(-K*(age-ao))*K*(age-ao)-Linf*exp(-K*(age-ao)),-Linf*exp(-K*(age-ao))*K*K)
  return(ret)
}


## hessian matrix - not a simple function!;

der2.vonb_nll = function(parm){
  Linf = exp(parm[1])
  K = exp(parm[2])
  ao = exp(parm[3])
  s2 = exp(parm[4])
  mu = vonb(Linf,K,ao,tdata$age)
  n = length(mu)
  resid = tdata$len-mu
  nll = n*log(2*pi*s2)/2 + sum(resid**2)/(2*s2)
  
  #second order derivatives 
  jac = der.vonb(Linf,K,ao,tdata$age)   
  dnll.dreg_parm = apply(-resid*jac,2,sum)/s2
  d2nll.dreg_parm = t(jac)%*%jac/s2
  
  #second order derivatives for s2
  t1 = apply(-resid*der.vonb.dLinf(Linf,K,ao,tdata$age)/s2,2,sum)
  t2 = apply(-resid*der.vonb.dK(Linf,K,ao,tdata$age)/s2,2,sum)
  t3 = apply(-resid*der.vonb.dao(Linf,K,ao,tdata$age)/s2,2,sum)
  d2nll.dreg_parm = d2nll.dreg_parm + cbind(t1,t2,t3);  
  
  d2nll.dreg_parm.ds2 = -dnll.dreg_parm/s2
  
  d2nll.d2s2 = -n/(2*s2*s2) + sum(resid**2)/(s2**3) 
  
  ret = matrix(NA,4,4)
  ret[1:3,1:3] = d2nll.dreg_parm
  ret[4,1:3]=d2nll.dreg_parm.ds2 
  ret[1:3,4]=d2nll.dreg_parm.ds2 
  ret[4,4]=d2nll.d2s2
  
  Ip = diag(exp(parm))
  
  ret = diag(der.vonb_nll(parm)) + Ip%*%ret%*%Ip 
  # note: this step because ret is deriv2 with respect to natural scale parms;
  # but we need deriv2s wrt to log(parm)
  
  rownames(ret) = c('Linf','K','ao','s2') 
  colnames(ret) = c('Linf','K','ao','s2') 
  
  return(ret)
}

## check
der2.vonb_nll(start.parm)

hessian(vonb_nll,start.parm)

#Newtons Algorithm to find the parameter values that minimize vonb_nll(parm);
parm=start.parm;
maxd = 1
iter=1

while(maxd>1e-6){
  invH = solve(der2.vonb_nll(parm));
  del = invH%*%der.vonb_nll(parm); # update in parm estimates
  parm = as.vector(parm - del); 
  fn = vonb_nll(parm)
  out = t(rbind(iter,fn,del)) 
  colnames(out) = c('iter','fn','del log(Linf)','del log(K)','del log(ao)','del log(s2)')
  print(out)
  print(exp(parm))
  maxd = max(abs(del))
  iter=iter+1
}

age.pred = seq(min(tdata$age),25,by=0.1)
mu.true = vonb(exp(parm.true[1]),exp(parm.true[2]),exp(parm.true[3]),age.pred); 
mu.fit = vonb(exp(parm[1]),exp(parm[2]),exp(parm[3]),age.pred);

plot(tdata$age,tdata$len,xlim=c(0,25),ylim = range(tdata$len,mu.true,mu.fit))
lines(age.pred,mu.true,lwd=2,col='blue')  
lines(age.pred,mu.fit,lwd=2,col='red')
legend("bottomright",bty='n',lwd=2,col=c('red','blue'),c('Fitted','True'))

##try different starting values;
parm=c(log(c(120,0.15,0.5,15)));
maxd = 1
iter=1
while(maxd>1e-6){
  invH = solve(der2.vonb_nll(parm));
  del = invH%*%der.vonb_nll(parm);
  parm = as.vector(parm - del); 
  fn = vonb_nll(parm)
  out = t(rbind(iter,fn,del)) 
  colnames(out) = c('iter','fn','del log(Linf)','del log(K)','del log(ao)','del log(s2)')
  print(out)
  print(exp(parm))
  maxd = max(abs(del))
  iter=iter+1
}

## better approach 
parm=c(log(c(120,0.15,1,15)));
maxd = 1
iter=1
while(maxd>1e-6){
  del = solve(der2.vonb_nll(parm),-der.vonb_nll(parm))
  parm = parm + del; 
  fn = vonb_nll(parm)
  out = c(iter,fn,del) 
  names(out) = c('iter','fn','del log(Linf)','del log(K)','del log(ao)','del log(s2)')
  print(out)
  print(exp(parm))
  maxd = max(abs(del))
  iter=iter+1
}

## better approach but poor starting values
parm=c(log(c(115,0.15,1,15)));
maxd = 1
iter=1
while(maxd>1e-6){
  del = solve(der2.vonb_nll(parm),-der.vonb_nll(parm))
  parm = parm + del; 
  fn = vonb_nll(parm)
  out = c(iter,fn,del) 
  names(out) = c('iter','fn','del log(Linf)','del log(K)','del log(ao)','del log(s2)')
  print(out)
  print(exp(parm))
  maxd = max(abs(del))
  iter=iter+1
}
## better approach but poor starting values; use a smaller newton step;
parm=c(log(c(115,0.15,1,15)));
maxd = 1
iter=1
while(maxd>1e-6){
  del = solve(der2.vonb_nll(parm),-der.vonb_nll(parm))
  parm = parm + 0.01*del; 
  fn = vonb_nll(parm)
  out = c(iter,fn,del) 
  names(out) = c('iter','fn','del log(Linf)','del log(K)','del log(ao)','del log(s2)')
  print(out)
  print(exp(parm))
  maxd = max(abs(del))
  iter=iter+1
}
opt.parm=parm

######################################So much easier with nlminb!!!
parm=c(log(c(115,0.15,1,15)));
fit1=nlminb(parm,vonb_nll)
fit1$par
opt.parm

# use gradient;
fit=nlminb(parm,vonb_nll,der.vonb_nll);
fit$par
opt.parm

