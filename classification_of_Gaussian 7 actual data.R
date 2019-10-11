#sample data
rm(list = ls())
getwd()
setwd("G:/")
library("ggplot2")
dat=readRDS("Adrenal-EMRE-dynamic.RDS")
ggplot(data = dat, aes(x = time, y = calc_mean, color = as.factor(malignancy) )) +       
  geom_line(aes(group = lesionID)) + geom_point()

dat0<-cbind.data.frame(mal=as.factor(dat$malignancy),
                       id=as.factor(dat$lesionID),
                       time=as.integer(dat$time),
                       mean=as.numeric(dat$calc_mean))
#reshape
library(reshape)
datz<-reshape(dat0, timevar="time", idvar=c("mal","id"), direction="wide")

##########Analysis 
require("fda.usc")
mdata=datz[3:dim(datz)[2]]
M=dim(mdata)[2]
head(mdata)

#dealing with missing data in arguments
#first count the missing values
table(is.na(mdata))
##replace with zoo
require(zoo)
#mdata<-na.approx(mdata, rule=2, na.rm=TRUE) #interpolate by columns
mdata<-t(na.approx(t(mdata), rule=2, na.rm=TRUE)) #interpolate by rows
colnames(mdata)<-paste("m",1:M,sep="_")
table(is.na(mdata))
head(mdata)
#interpolating missing values is not possible with zoo alone
##replace with lesion average
miss<-which(is.na(mdata), arr.ind=TRUE)
mdata[miss]<-rowMeans(mdata, na.rm=TRUE)[miss[,1]]
table(is.na(mdata))
head(mdata)

##true analysis begins here
mlearn=fdata(mdata)
glearn=datz$mal
dataf<-data.frame(glearn) 
dat1=list("df"=dataf,"x"=mlearn)
a1<-classif.glm(glearn~x, data=dat1)
newdat<-list("x"=mlearn) 
p1<-as.factor(predict.classif(a1,newdat) )
acc=function(x, y){sum(x==y)/length(x)}
err.glm=acc(p1, glearn)
err.glm

out.glm=classif.glm(glearn~x, data = dat1)
out.np=classif.np(glearn,mlearn)
out.knn=classif.knn(glearn,mlearn,knn=c(3,5,7))
out.kernel=classif.kernel(glearn,mlearn)

res.glm=summary.classif(out.glm)
res.np=summary.classif(out.np)
res.knn=summary.classif(out.knn)
res.kernel=summary.classif(out.kernel)

res.glm$prob.classification
res.np$prob.classification
res.knn$prob.classification
res.kernel$prob.classification

#Machine learning models
require(caret)
#neural network
nn<-train(
  x=mdata,
  y=datz$mal,
  method="nnet",
  trace=FALSE)
nn
confusionMatrix(nn)
#random forests
rf<-train(
  x=mdata,
  y=datz$mal,
  method="ranger")
rf
confusionMatrix(rf)
