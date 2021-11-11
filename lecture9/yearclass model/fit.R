setwd("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week8\\3LNO_ex")

source("C:\\home\\CADIGAN\\GradProgram\\2018\\F6004\\week2\\spay.txt")


#install.packages('Matrix') 
#install.packages('latticeExtra')
#install.packages('corrplot')     
#install.packages('reshape2')             
#install.packages('lme4')    
#install.packages('cAIC4')    
#install.packages('xtable')   
#install.packages('ggplot2') 
  
require(lattice)
library(latticeExtra)
require(reshape2)  
library(lme4)          
library(cAIC4)  
library(xtable)  
library(ggplot2)  
library(psych)    
library(stargazer)

sdat = read.table(file='F3LNO.dat',header=TRUE)

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

colv = c("blue","cyan","aquamarine4","green","darkgoldenrod4","red")

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Index"), cex=1)
stripttl <- list(cex=0.75,font=1) 
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
my.key <- list(x = .96, y = .96, corner = c(1, 1),
            size = 2,between=0.5,
            text = list(as.character(0:5),cex=0.6),
            lines = list(type = c("l"), col = colv,
                         lwd = 2, lty = 1,span=0.25))
 
P1 = xyplot(index~Year,groups=Age,cex=1,
data=sdat, type=c("l"),lty=1,xlab=xttl,lwd=2,
ylab=yttl, strip.text=stripttl, las=1,
key=my.key,col=colv,par.settings = my.padding,
par.strip.text = list(cex=0.6),
) 
print(P1)

gname=c('RV.jpeg')
jpeg(file=gname,width=3,height=3,units='in',res=300)
par(mar=c(2.5,3,0,0),mgp=c(2,0.75,0))
print(P1) 
dev.off()

sdat$log.index.dev=resid(glm(log.index ~ factor(Age) -1,data=sdat))

my.key <- list(corner = c(1,0),
            size = 2,between=0.5,
            text = list(as.character(0:5),cex=0.6),
            lines = list(type = c("l"), col = colv,
                         lwd = 2, lty = 1,span=0.25))

xttl <- list(c("Cohort"), cex=1)
yttl <- list(c("Log Index Deviation"), cex=1)
P2 = xyplot(log.index.dev~YC,groups=Age,cex=1,
data=sdat, type=c("l"),lty=1,xlab=xttl,lwd=2,
ylab=yttl, strip.text=stripttl, las=1,
key=my.key,col=colv,par.settings = my.padding,
par.strip.text = list(cex=0.6),
) 
print(P2)

gname=c('log_RV_dev.jpeg')
jpeg(file=gname,width=3,height=3,units='in',res=300)
par(mar=c(2.5,3,0,0),mgp=c(2,0.75,0))
print(P2) 
dev.off() 


t1=dcast(sdat[,2:4],Year~Age)

scale=1
 
sp.dat <- t(spay(t(t1[,2:7]))) 
syear <- rep(t1[,1],length(0:6))
sage <- rep(0:6,each=length(t1[,1]))
vspay <- as.vector(unlist(sp.dat))

cohort.ps=c(0.02,0.01,0.03,0.01)
y.adj = 0.2
ylim=range(sage);
ylim[2]=ylim[2]+y.adj
xlim=range(syear);
xlim[2]=xlim[2]+0.2

gname <- 'spay.jpeg'
jpeg(file=gname,width=3,height=3,units='in',res=300)

par(mar=c(3,3,1.2,0.2),cx.axis=0.7)

ind = !is.na(vspay)
bp(syear[ind],sage[ind],vspay[ind],cohort.ps=cohort.ps, ylim=ylim, xlim=xlim, 
     las=1, xlab='', ylab='', main="", scale=scale)
points(syear[ind],sage[ind],cex=sqrt(abs(vspay[ind]))*scale, col='black', pch=1)
mtext(side=3,outer=F,line=.1,c('3NO Cod Fall RV (ages 0-5) SPAY'),cex=0.8) 
mtext(side=2,outer=F,line=2,'Age')  
mtext(side=1,outer=F,line=2,'Year')
  
dev.off()

##################   Fixed Effects Model ##############


sdat$fYC = factor(sdat$YC) 
sdat$fAge = factor(sdat$Age)
glmfit = glm(log.index ~ fYC + fAge -1,data=sdat)

sdat$pred = predict(glmfit)
sdat$resid = resid(glmfit)  
sdat$sresid = rstandard(glmfit)

stargazer(glmfit, type = "html",  out="out.doc",align = TRUE,   
 title="Cohort Strength, Fall Survey 3LNO cod", single.row=TRUE)
# ci=TRUE, ci.level=0.95)
    
out = data.frame(log_est=coef(glmfit),ci=confint(glmfit))
rec.dev = subset(out,substr(rownames(out),1,3)=='fYC')  
rnames = rownames(rec.dev)
rec.dev$YC = as.numeric(substring(rnames,4,7))
rec.dev$est = exp(rec.dev$log_est)   
mean.est = mean(rec.dev$est)  
rec.dev$est = rec.dev$est/mean.est   
rec.dev$L95 = exp(rec.dev$ci.2.5..)/mean.est 
rec.dev$U95 = exp(rec.dev$ci.97.5..)/mean.est


ylim = range(rec.dev[,5:7],na.rm=T)
ylim[2]=4

gname = "YC.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)
  par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(est~YC, data=rec.dev,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Cohort") 
  mtext(side=2,line=2.3,"Relative Cohort Strength")
  
  sx = rec.dev$YC
  low = rec.dev$L95
  high = rec.dev$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),col = grey(0.7),border = NA)
  
  
  vy = unique(sdat$YC)
  vy = vy[2:(length(vy)-1)]
  abline(v=vy,lty=2,col="darkgoldenrod1")
  abline(h=1,lty=2)
  lines(rec.dev$YC,rec.dev$est,lwd=2) 

dev.off()


fit=sdat
fit$survey="Fall 3LNO"
fit$resid = sdat$sresid

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
lines(unique(x),mres,lty=1,lwd=1,col='red')
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

colv

sdat$colv = c("a","b","c","d","e")[sdat$YC-5*floor(sdat$YC/5)+1]
sdat$block=1
sdat$block[sdat$YC>2003]=2 
sdat$block[sdat$YC>2008]=3

mytheme <- theme(panel.background = element_blank(),
                 axis.line=element_line(color="gray"),
                 panel.grid.major = element_line(color="lightgray", linetype = 3,size=.2),
                 panel.grid.minor = element_line(color="lightgray", linetype = 3,size=.2),
                 panel.border = element_rect(colour = "black",fill=NA, size=.1),
                 axis.text.x=element_text(angle=0,hjust=0,vjust=0,size=8),
                 plot.title=element_text(hjust=0.5),
                 legend.text=element_text(size=6),
               #  legend.title = element_text(size=4),
                 legend.title=element_blank(),
                 legend.key.width = unit(0.7,"line"),
                 legend.key=element_rect(colour = "white", fill = "white"),  
                 legend.position="bottom")
                 
plot2 <-ggplot(data=sdat,aes(Year,pred,fill=factor(Age),color=factor(colv)))+
  geom_point(shape=21,size=1,aes(y=log.index),alpha=0.4,show.legend=F)+
  geom_line(aes(y=log.index,group=YC,color=colv),alpha=0.6,show.legend=F)+
  geom_line(aes(group=YC),size=0.7,color='black')+
  
  mytheme + 
#  scale_fill_brewer("YC",type="seq",palette="GrRdBl")+ 
 # scale_fill_manual(values=colv)+
  scale_color_manual(values=c("red","blue","green","orange","aquamarine4"))+
  theme(legend.margin=margin(t = -0.6, unit='cm'),legend.position = c(0.8, .2),
  legend.justification = c("right", "top"),legend.key.height=unit(0.5,"line"),
  strip.text.x = element_text(size=7,colour = "black", angle = 0, margin = margin(0.5,0,0.5,0, "mm")))+
  facet_wrap(~block,ncol=1,scales = "free_x")+  
 # guides(fill=FALSE)+       
guides(color=guide_legend(ncol=4))+
  labs(x = "Year",y='Log Survey Index',color="Ages",fill="Ages")            

print(plot2)     
                 
gname=c('CC.jpeg')
jpeg(file=gname,width=3.5,height=5,units='in',res=300)       
  print(plot2)
dev.off()


sdat$colr = 'deepskyblue'
sdat$colr[sdat$resid>0]='firebrick2'

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year , sdat, 
       ylab='Age',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 3*sqrt(abs(sdat$sresid)/pi), fill.color = sdat$colr,
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

################# mixed effects model ############################################;

sdat$fYear = factor(sdat$Year)                                                                                                    
mixfit <- lmer(log.index ~ fYC + fAge -1 + (1|fYear),REML = FALSE, data=sdat)


sdat$mix.pred = predict(mixfit)     
sdat$mix.resid = residuals(mixfit,type="response")
sdat$mix.sresid = residuals(mixfit,type="pearson") 

stargazer(mixfit, type = "html",  out="out1.doc",align = TRUE,   
 title="Cohort Strength, Fall Survey 3LNO cod", single.row=TRUE)
# ci=TRUE, ci.level=0.95)

stargazer(ranef(mixfit)$fYear, type = "html",  out="out2.doc",align = TRUE,   
 title="Year Effects, Fall Survey 3LNO cod", single.row=TRUE,summary=FALSE)
 
hist(ranef(mixfit)$fYear[,1],breaks=10)
fYear = ranef(mixfit)$fYear[,1]
std.fyear = data.frame(VarCorr(mixfit))[1,5]
s.Year = fYear/std.fyear
qqnorm(s.Year)
abline(a=0,b=1)

year = as.numeric(rownames(ranef(mixfit)$fYear))

par(mfrow=c(2,1),mar=c(3,3,0.5,0.5),mgp=c(2,0.5,0))
plot(year,fYear,xlab='Year',ylab='Year Effect',type='b')
abline(h=0,lty=2)

acf(fYear,plot=TRUE)

pname = names(fixef(mixfit))
 
fcohort = fixef(mixfit)[substr(pname,1,3)=='fYC']
fcohort  = fcohort-mean(fcohort)
cohort = as.numeric(gsub('fYC','',names(fcohort))) 

par(mfrow=c(2,1),mar=c(3,3,1,0.5),mgp=c(2,0.5,0))
plot(cohort,fcohort,xlab='Cohort',ylab='Rel. Cohort Size',type='b')
abline(h=0,lty=2)

acf(fcohort,plot=TRUE)

fcohort.resid = resid(glm(fcohort~cohort))

par(mfrow=c(2,1),mar=c(3,3,1,0.5),mgp=c(2,0.5,0))
plot(cohort,fcohort.resid,xlab='Cohort',ylab='Rel. Cohort Size',type='b')
abline(h=0,lty=2)

acf(fcohort.resid,plot=TRUE)


ci=confint(mixfit,parm=names(fixef(mixfit)),method = "Wald")    
out = data.frame(log_est=fixef(mixfit),ci=ci)
mix.rec.dev = subset(out,substr(rownames(out),1,3)=='fYC')  
rnames = rownames(mix.rec.dev)
mix.rec.dev$YC = as.numeric(substring(rnames,4,7))
mix.rec.dev$est = exp(mix.rec.dev$log_est)   
mean.est = mean(mix.rec.dev$est)  
mix.rec.dev$est = mix.rec.dev$est/mean.est   
mix.rec.dev$L95 = exp(mix.rec.dev$ci.2.5..)/mean.est 
mix.rec.dev$U95 = exp(mix.rec.dev$ci.97.5..)/mean.est


ylim = range(rec.dev[,5:7],mix.rec.dev[,5:7],na.rm=T)
ylim[2]=4

gname = "YC_mixed.jpeg"
jpeg(file=gname,width=3,height=3,units='in',res=300)
  par(mar=c(3,3.3,0.2,1),mgp=c(2,0.7,0))
  plot(est~YC, data=rec.dev,ylab = '', xlab='',las=1,type='l',lwd=2,ylim=ylim)
  mtext(side=1,line=1.7,"Cohort") 
  mtext(side=2,line=2.3,"Relative Cohort Strength")
  
  sx = rec.dev$YC
  low = rec.dev$L95
  high = rec.dev$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),col = grey(0.7),border = NA)
  sx = mix.rec.dev$YC
  low = mix.rec.dev$L95
  high = mix.rec.dev$U95
  polygon(c(sx,rev(sx)),c(low,rev(high)),col = rgb(0,0,1, alpha=0.2, maxColorValue = 1),border = NA)
  
  vy = unique(sdat$YC)
  vy = vy[2:(length(vy)-1)]
  abline(v=vy,lty=2,col="darkgoldenrod1")
  abline(h=1,lty=2)
  lines(rec.dev$YC,rec.dev$est,lwd=2) 
  lines(mix.rec.dev$YC,mix.rec.dev$est,lwd=2,col='blue')
  
  legend("topleft",col=c('black','blue'),lty=1,lwd=2,legend=c('glm','lmer'),bty='n')

dev.off()

fit=sdat
fit$survey="Fall 3LNO"
fit$resid = sdat$mix.sresid

gname <- "sres_mix.jpeg"
jpeg(file=gname,width=4,height=5,units='in',res=300)
par(mfcol=c(4,1),oma=c(0,3,1,1),mar=c(3,2,0,2),las=1)

ind <- !is.na(fit$pred)
x <- fit$Year[ind]
y <- fit$resid[ind]

plot(x,y,xlab="",ylab="",type='n')
text(x,y,fit$Age[ind])
abline(h=0,lty=1)
mres <- tapply(y,x,"mean")
lines(unique(x),mres,lty=1,lwd=1,col='red')
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

              
plot2 <-ggplot(data=sdat,aes(Year,mix.pred,fill=factor(Age),color=factor(colv)))+
  geom_point(shape=21,size=1,aes(y=log.index),alpha=0.4,show.legend=F)+
  geom_line(aes(y=log.index,group=YC,color=colv),alpha=0.6,show.legend=F)+
  geom_line(aes(group=YC),size=0.7,color='black')+
  
  mytheme + 
#  scale_fill_brewer("YC",type="seq",palette="GrRdBl")+ 
 # scale_fill_manual(values=colv)+
  scale_color_manual(values=c("red","blue","green","orange","aquamarine4"))+
  theme(legend.margin=margin(t = -0.6, unit='cm'),legend.position = c(0.8, .2),
  legend.justification = c("right", "top"),legend.key.height=unit(0.5,"line"),
  strip.text.x = element_text(size=7,colour = "black", angle = 0, margin = margin(0.5,0,0.5,0, "mm")))+
  facet_wrap(~block,ncol=1,scales = "free_x")+  
 # guides(fill=FALSE)+       
guides(color=guide_legend(ncol=4))+
  labs(x = "Year",y='Log Survey Index',color="Ages",fill="Ages")            

print(plot2)     
                 
gname=c('CC_mix.jpeg')
jpeg(file=gname,width=3.5,height=5,units='in',res=300)       
  print(plot2)
dev.off()


sdat$colr = 'deepskyblue'
sdat$colr[sdat$mix.resid>0]='firebrick2'

jDarkGray <- 'grey20'
jPch <- 21
resid.plot1 = xyplot(factor(Age)~Year , sdat, 
       ylab='Age',
       scales = list(y = 0:5,x=c(2000,2005,2010,2015)),
       cex = 3*sqrt(abs(sdat$mix.sresid)/pi), fill.color = sdat$colr,
       col = jDarkGray,
       par.settings = my.padding,
       par.strip.text = list(cex=0.5),
       panel = function(x, y, ..., cex, fill.color, subscripts) { 
         panel.abline(h=1:6, v=c(2000,2005,2010,2015), col.line=grey(0.9))
         panel.xyplot(x, y, cex = cex[subscripts],
                      pch = jPch, fill = fill.color[subscripts], ...)
         })
         
print(resid.plot1)
                  
gname = c("resid_matrix_bubbles_mix.jpeg")
jpeg(file=gname,width=4,height=3,units='in',res=300)
  print(resid.plot1)
dev.off()








