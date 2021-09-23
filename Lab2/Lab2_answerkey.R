#set working directory 
setwd("C:/Users/aperreau/Google Drive/Fall 2021/FISH6004_Lab/Lab2/")

#read in data
Data<- read.csv(file = "ALK_sample.csv")
Data_catch <- read.csv(file = "catch_num.csv")

head(Data)
head(Data_catch)

#Part I
#calculate the proportion at age
n_ages = ncol(Data)-2
n_total = Data[,ncol(Data)]
n_total_mat = matrix(unlist(rep(n_total, n_ages)), ncol=n_ages)
proportion = Data[,2:(n_ages+1)]/n_total_mat

#check if this looks right
print(proportion)

#rows should sum to 1
print(apply(proportion,1,FUN = sum))

#plot age at length
par(mfcol=c(1,1),mar=c(3,3,2,2),las=1,mgp=c(2,1,0),oma=c(1,1,1,1))
image(as.matrix(proportion), xaxt="n", yaxt= "n", xlab="Length", ylab="Age")
axis(1,at = seq(0, 1, length.out = 20), labels = 1:20)
axis(2,at = seq(0, 1, length.out = 25), labels = 1:25)

#Calculate catch_at_age_length and catch at age
Data_catch_mat <- matrix(unlist(rep(Data_catch[,2],n_ages)), ncol = n_ages)
catch_at_age_length <- Data_catch_mat*proportion
catch_at_age <- apply(catch_at_age_length,2,FUN = sum)
names(catch_at_age) <- 1:25

catch_at_length <- apply(catch_at_age_length,1,FUN = sum)
names(catch_at_length) <- 1:20
par(mfcol=c(1,2),mar=c(3,3,2,2),las=1,mgp=c(2,1,0),oma=c(1,1,1,1))
barplot(catch_at_age, xlab="Age", ylab="Catch")
barplot(catch_at_length, xlab="Length", ylab="Catch")

#####################
#Part II: Use catch at age data from 3Ps cod to generate a spay plot

#Read in the data.
data <- read.csv("3pscod_C@A.csv")

#Plot the age composition data
n.year = nrow(data)
n.age = ncol(data)-1
Cm = data[,2:ncol(data)]/1000 
Bs = t(Cm) 
colnames(Bs) = data[,1]
rownames(Bs) = as.character(3:14)

col1 <- colorRampPalette(rev(c("red","green","yellow","cyan", "blue")))  
barplot(Bs, col=col1(nrow(Bs)) , border="white", space=0.04, cex.axis=0.5,cex.names=0.5,
        xlab="Year",ylab='',main='') 
mtext(side=2,line=2,"Proportion at age 1-10+",las=0,outer=T)
mtext(side=3,line=0,"Numbers at age (000's)",las=0)

Bs1 = Bs%*%diag(1/apply(Bs,2,sum)) 
colnames(Bs1) = data[,1]
barplot(Bs1, col=col1(nrow(Bs1)) , border="white", space=0.04, cex.axis=0.5,cex.names=0.5,
        xlab="Year",ylab='',main='')
legend("topright", rownames(Bs), fill = col1(nrow(Bs)), bty = "n",
       xjust=1,inset=c(-0.09,0),xpd=NA,cex=0.7,x.intersp=0.1) 
mtext(side=3,line=0,"Proportion at age",las=0)

#generate a 'spay' plot. 
# define functions from Noel

#this will generate the spay plot
bp<-function(x,y,v, scale=3,cohort.id=T,cohort.ps, ...){
  plot(x,y,cex=sqrt(abs(v))*scale, col=ifelse(v<0,'black','grey'), pch=ifelse(v<0,16,16), ...)
  points(x[v>0],y[v>0],cex=sqrt(v[v>0])*scale, col='grey', pch=1, ...)
  
  if(cohort.id){
    lx <- max(x)             
    yval <- unique(y) 
    xval <- unique(x) 
    yval1 <- yval[1:(length(yval)-1)]
    surv.cohort<- lx-yval1
    surv.cohort1 <- substring(as.character(surv.cohort),3,4)
    pos1 <- lx + cohort.ps[1]*diff(range(xval))
    pos2 <- yval1 + cohort.ps[2]*diff(range(yval1))
    text(pos1,pos2,surv.cohort1,cex=0.5,srt=45)
    
    xval1 <- xval[seq(2,length(xval),by=2)]
    ly <- max(y)
    pos1 <- ly + cohort.ps[3]*diff(range(yval))
    pos2 <- xval + cohort.ps[4]*diff(range(xval))
    old.cohort<- xval-ly
    old.cohort1 <- substring(as.character(old.cohort),3,4)
    text(pos2,pos1,old.cohort1,cex=0.5,srt=45)
  }  
}

#this calculates standardized proportions at age per year
spay <- function(tdat){
  stdat <- apply(tdat,2,sum,na.rm=T)
  ny <- length(stdat)
  na <- nrow(tdat)
  mstdat <- matrix(stdat,nrow=na,ncol=ny,byrow=T)
  ptdat <- tdat/mstdat
  ptdat.dev <- ptdat - matrix(apply(ptdat,1,mean,na.rm=T),nrow=na,ncol=ny,byrow=F)
  ptdat.std <- sqrt(apply(ptdat.dev^2,1,mean,na.rm=T))
  ztdat <- ptdat.dev/matrix(ptdat.std,nrow=na,ncol=ny,byrow=F)
  return(ztdat)}

scale=1.3

#
sdat <- t(spay(t(Cm))) 
syear <- rep(data[,1],n.age)
sage <- rep(3:14,each=n.year)
vspay <- as.vector(unlist(sdat))

cohort.ps=c(0.02,0.01,0.03,0.01)
y.adj = 0.2
ylim=range(sage);
ylim[2]=ylim[2]+y.adj
xlim=range(syear);
xlim[2]=xlim[2]+0.2

par(mfcol=c(1,1),mar=c(3,3,2,2),las=1,mgp=c(2,1,0),oma=c(1,1,1,1))
bp(syear,sage,vspay,cohort.ps=cohort.ps, ylim=ylim, xlim=xlim, 
   las=1, xlab='', ylab='', main="", scale=scale)
points(syear,sage,cex=sqrt(abs(vspay))*scale, col='black', pch=1)
mtext(side=3,outer=F,line=.1,c('3Ps Cod C@A SPAY'),cex=1) 
mtext(side=2,outer=F,line=2,'Age',cex=1)  
mtext(side=1,outer=F,line=2,'Year',cex=1)


#replicate the spay plot using ggplot 
library(reshape2)
library(tidyverse)

temp<- Cm
rownames(temp) = c(1959:2016)
colnames(temp) = c(3:14)
temp_vec<-melt(as.matrix(temp))
colnames(temp_vec) = c("Year","Age","catch")
temp_vec$YC <- temp_vec$Year-temp_vec$Age 

sdat.spay <- temp_vec%>%group_by(Year)%>%mutate(prop.index=catch/sum(catch,na.rm=T))%>%
  group_by(Age)%>%mutate(std.prop.index=scale(prop.index,scale=T))%>%ungroup()%>%
  mutate(col=ifelse(std.prop.index<0,"-Ve","+Ve"),
         text=ifelse(Age==max(Age)|Year==max(Year),YC,NA),
         text1=substr(text,3,4))

ggplot(sdat.spay,aes(x=Year,y=Age,size=std.prop.index,fill=col))+
  geom_point(shape=21,alpha=0.8,show.legend=F)+
  scale_fill_manual(values=c("black","grey"))+
  scale_size_area()+
  geom_text(aes(label=text1),hjust=-1,angle=45,color="black", size=2.2)+
  scale_x_continuous(breaks = seq(1960,2018,8))+
  scale_y_continuous(breaks = seq(0,15,1))+
  labs(x="Year",y="Age",title="3Ps Cod C@A SPAY")+
  theme(text = element_text(size=5))+
  theme_bw()

