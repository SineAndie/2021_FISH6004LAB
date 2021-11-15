require(reshape2) 
require(lattice)
library(latticeExtra)

dir.name = getwd()

source("spay.txt")
 

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

cnames = c('Year',paste('Age',0:12,sep=""))
FRV.matrix = read.table(file='data\\FRV.dat',header=FALSE,col.names=cnames)
SRV.matrix = read.table(file='data\\SRV.dat',header=FALSE,col.names=cnames) 
cnames = c('Year',paste('Age',0:14,sep=""))       
JRV.matrix = read.table(file='data\\JRV.dat',header=FALSE,col.names=cnames)

cnames = c('Year',paste('Age',2:12,sep=""))
catch = read.table(file='data\\catch.dat',header=FALSE,col.names=cnames) 

cnames = c('Year',paste('Age',1:12,sep=""))
mat = read.table(file='data\\mat.txt',header=FALSE,col.names=cnames)

cnames = c('Year',paste('Age',3:12,sep=""))
wt = read.table(file='data\\wt.txt',header=FALSE,col.names=cnames) 

## subset for assessment

FRV.matrix = FRV.matrix[,c(1,4:12)]
SRV.matrix = SRV.matrix[,c(1,4:12)] 
JRV.matrix = JRV.matrix[,c(1,4:12)]

mat = mat[,3:13]
wt[,1] = 0
colnames(wt)=colnames(mat)

year = catch[,1]
catch = catch[,2:12]/1000

# make vector data frames

vec_func = function(mdat){
  vdat = melt(mdat,id=c("Year"))
  vdat$Age = as.numeric(substr(vdat$variable,4,5))
  vdat$variable=NULL
  vdat$index=vdat$value
  vdat$value=NULL
  return(vdat)
}

## make a spay plot

plot.spay = function(x,year,age,yname,gname,x.size,y.size,scale){

  sp.dat <- t(spay(t(x))) 
  syear <- rep(year,length(2:12))
  sage <- rep(2:12,each=length(x[,1]))
  vspay <- as.vector(unlist(sp.dat))

  cohort.ps=c(0.02,0.01,0.03,0.01)
  y.adj = 0.2
  ylim=range(sage);
  ylim[2]=ylim[2]+y.adj
  xlim=range(syear);
  xlim[2]=xlim[2]+0.2

  #gname <- 'data\\catch_spay.jpeg'
  jpeg(file=gname,width=x.size,height=y.size,units='in',res=300)

  par(mar=c(3,3,1.2,0.2),cex.axis=0.7)

  ind = !is.na(vspay)
  bp(syear[ind],sage[ind],vspay[ind],cohort.ps=cohort.ps, ylim=ylim, xlim=xlim, 
     las=1, xlab='', ylab='', main="", scale=scale)
  points(syear[ind],sage[ind],cex=sqrt(abs(vspay[ind]))*scale, col='black', pch=1)
  mtext(side=3,outer=F,line=.1,yname,cex=0.8) 
  mtext(side=2,outer=F,line=2,'Age')  
  mtext(side=1,outer=F,line=2,'Year')
  
  dev.off() 
}


FRV.vec = vec_func(FRV.matrix)
FRV.vec$survey = 'Fall'     
SRV.vec = vec_func(SRV.matrix)
SRV.vec$survey = 'Spring'    
JRV.vec = vec_func(JRV.matrix)
JRV.vec$survey = 'Juvenile' 

indices = rbind(FRV.vec,SRV.vec,JRV.vec)

temp = cbind(year,catch)
colnames(temp) = c('Year',colnames(catch))
catch.vec = vec_func(temp) 

temp = cbind(c(year,max(year)+1),wt)
colnames(temp) = c('Year',colnames(wt))
wt.vec = vec_func(temp)    

temp = cbind(c(year,max(year)+1),mat)
colnames(temp) = c('Year',colnames(mat))
mat.vec = vec_func(temp)

##plot inputs

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Catch (millions)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
dat = catch.vec 
nr = ceiling(length(2:12)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\catch_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()

plot.spay(catch,year,age,'3NO Cod Catch (ages 0-12) SPAY','data\\catch_spay.jpeg',6,4,1)

ind = year>=1994
plot.spay(catch[ind,],year[ind],age,'3NO Cod Catch (ages 0-12) SPAY','data\\catch_spay_recent.jpeg',4,4,1.5)


## plot Spring Survey

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Spring RV Index (mnpt)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
dat = SRV.vec 
nr = ceiling(length(2:10)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\SRV_ts.jpeg"
jpeg(file=gname,width=4,height=3,units='in',res=300)       
  print(P1)
dev.off()

plot.spay(SRV.matrix[,2:10],SRV.matrix[,1],2:10,'3NO Cod Spring RV SPAY','data\\SRV_spay.jpeg',4,4,1)

## plot Fall Survey

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Fall RV Index (mnpt)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
dat = FRV.vec 
nr = ceiling(length(2:10)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\FRV_ts.jpeg"
jpeg(file=gname,width=4,height=3,units='in',res=300)       
  print(P1)
dev.off()

plot.spay(FRV.matrix[,2:10],FRV.matrix[,1],2:10,
'3NO Cod Fall RV SPAY','data\\FRV_spay.jpeg',4,4,1)


## plot Fall Survey


xttl <- list(c("Year"), cex=1)
yttl <- list(c("Juvenile RV Index (mnpt)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=1,relation="free",tck=0.3,y=list(rot=0))
dat = JRV.vec 
nr = ceiling(length(2:10)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\JRV_ts.jpeg"
jpeg(file=gname,width=4,height=3,units='in',res=300)       
  print(P1)
dev.off()

plot.spay(JRV.matrix[,2:10],JRV.matrix[,1],2:10,
'3NO Cod Juvenile RV SPAY','data\\JRV_spay.jpeg',4,4,2)

#Plot mats

xttl <- list(c("Year"), cex=1)
yttl <- list(c("% Mature"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.6,relation="free",tck=0.3,y=list(rot=0))
dat = mat.vec
dat$index = round(100*dat$index,digits=1) 
nr = ceiling(length(2:12)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\mat_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()  

#Plot weights

xttl <- list(c("Year"), cex=1)
yttl <- list(c("Body Weight (kg)"), cex=1)
stripttl <- list(cex=0.5,font=1)
ax <- list(cex=0.6,relation="free",tck=0.3,y=list(rot=0))
dat = wt.vec
#dat$index = round(100*dat$index,digits=1) 
nr = ceiling(length(2:12)/3)

P1 = xyplot(index~Year|as.factor(Age),
data=dat, type=c("l"),lwd=2,xlab=xttl, scales=ax,
ylab=yttl, strip.text=stripttl, las=1,col=c('black'),
as.table=TRUE, layout=c(3,nr,1),
par.settings = my.padding,
par.strip.text = list(cex=0.6),
)
gname = "data\\wt_ts.jpeg"
jpeg(file=gname,width=4,height=4,units='in',res=300)       
  print(P1)
dev.off()


save(catch,mat,wt,indices,catch.vec,mat.vec,wt.vec,file='data//3LNO.RData')

