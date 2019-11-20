rm(list=ls())
setwd("E:/2. Other research/1. Classification of Gaussian Process/actual data")

library(ggplot2)
library(zoo)
library(fda.usc)
library(tidyr)
library(dplyr)

#load data set
dat=readRDS("Adrenal-EMRE-dynamic.RDS")
head(dat)

#separate out the working data and rename variables
dat0<-cbind.data.frame(mal=as.factor(dat$malignancy),
                       id=as.factor(dat$lesionID),
                       time=as.integer(dat$time),
                       mean=as.numeric(dat$calc_mean))
head(dat0)

#plot the data set
ggplot(data=dat0,aes(x=time,y=mean,color=mal))+
  geom_point()+
  geom_line(aes(group=id), alpha = 0.2)+
  ggtitle("Scan results over time")+
  xlab("Time of scan")+
  ylab("Mean Index")+
  labs(color="Malignancy")+
  theme_bw()


#reshape
dat1<-spread(dat0,time,mean)
head(dat1)

#Dealing with the missing data
#choose only the data matrix
#dat2<-dat1[, grepl("mean",names(dat1))]

#replace missing with mean
dat2<-dat1
dat2_miss<-which(is.na(dat2),arr.ind=TRUE)
dat2[dat2_miss]<-rowMeans(dat2[,3:dim(dat2)[2]],na.rm=TRUE)[dat2_miss[,1]]
head(dat2)

#replace missing using zoo constant method
dat2<-dat1
for(i in 1:dim(dat2)[1]){
  temp<-c(dat2[i,3:dim(dat2)[2]])
  temp<-na.approx(temp,rule=2,na.rm=FALSE)
  dat2[i,3:dim(dat2)[2]]<-temp
}
head(dat2)

#replace missing using zoo constant method
dat2<-dat1
for(i in 1:dim(dat2)[1]){
  temp<-c(dat2[i,3:dim(dat2)[2]])
  temp<-na.approx(temp,rule=2,na.rm=FALSE)
  dat2[i,3:dim(dat2)[2]]<-temp
}
head(dat2)

#make the new data long format and compare
dat3<-gather(dat2,"time","mean",3:dim(dat2)[2])
dat4<-merge(dat0,dat3,by=c("id","time","mal"),all=TRUE)
dat4$time<-as.numeric(dat4$time)
head(dat4)

ggplot(data=dat4)+
  geom_point(aes(x=time,y=mean.x,color=mal))+
  geom_line(aes(x=time,y=mean.y,group=id),alpha=.1)+
  ggtitle("Scan results over time")+
  xlab("Time of scan")+
  ylab("Mean Index")+
  labs(color="Malignancy")+
  theme_classic()
