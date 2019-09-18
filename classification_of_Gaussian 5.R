rm(list=ls())
require(MASS)
library(ggplot2)
#install.packages("fda.usc")
#require(fda.usc)

#number of screening
T=5
#T=10
t=seq(0,2*(T-1),2)
#t=seq(1,T,1)
#number of patients
n=1500

#mean when not malignant
#mu0=.5*sin(.5*t)
mu0=10+log(1+1.5*t)
#mean when malignant
#mu1=.1+.6*sin(.5*t+.5)
mu1=10.1+log(1+2*t)

#plot means only
g1<-ggplot()+
  geom_line(aes(x=t,y=mu0),colour="red")+
  geom_line(aes(x=t,y=mu1),colour="blue")
g1

#probability of cancer for each patient
pc=runif(n)
#create binary response variable
malignant1<-as.numeric(pc>.5)
malignant<-c()
malignant[malignant1==1]<-"malignant"
malignant[malignant1==0]<-"non-malignant"

#create random noise for each patient
#the following function generates positive definite matrix
Posdef <- function (n, ev = runif(n, 0, 1)) 
{Z<-matrix(ncol=n, rnorm(n^2)); decomp<-qr(Z); Q<-qr.Q(decomp);
R<-qr.R(decomp); d<-diag(R); ph<-d/abs(d); O<-Q%*%diag(ph);
Z<-t(O)%*%diag(ev)%*%O; return(Z)}

Sigma=Posdef(T)

#create gaussian process for each patient
xs<-matrix(c(0),nrow=n,ncol=T)
for(i in 1:n){
  if(malignant[[i]]=="non-malignant"){
    xs[i,]=mvrnorm(1,mu=mu0,Sigma=Sigma)
  }
  if(malignant[[i]]=="malignant"){
    xs[i,]=mvrnorm(1,mu=mu1,Sigma=Sigma)
  }
}

#create id and data frame
t1=data.frame(rep(t,n))
xs1=data.frame(c(xs))
m1=data.frame(c(t(matrix(c(rep(malignant,T)),nrow=n,ncol=T,byrow=F))))
id1=data.frame(c(t(matrix(rep(rep(1:n),T),nrow=n,ncol=T,byrow=F))))
data=cbind.data.frame(m1,xs1,t1,id1)
colnames(data)<-c("m1","xs1","t1","id")

#see the data
head(data,10)

#plot data
g1+
  geom_point(aes(x=data$t1,y=data$xs1,colour=data$m1),alpha=.2)+
  labs(title="Test results", x="Test (time)", y="Malignancy")+
  scale_colour_manual(name="Malignancy",
                      values=c("malignant"="blue","non-malignant"="red"))

#Analysis
data1=cbind.data.frame(xs,t,malignant1,seq(1,n,1))
names(data1)[1:T] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep=""), x="X", y=T))
col_temp1<-T+1;col_temp2<-dim(data1)[2]
names(data1)[col_temp1:col_temp2]<-c("t","m","id")
head(data1,4)


require("fda.usc")
mdata= xs                          #matrix(xs[,1],byrow=T,ncol=length(t) )
mlearn=fdata(mdata,argvals=as.integer(t))
glearn=as.factor(malignant)
dataf<-data.frame(glearn) 
dat0=list("df"=dataf,"x"=mlearn)
a1<-classif.glm(glearn~x, data = dat0)
newdat<-list("x"=mlearn) 
p1<-as.factor(predict.classif(a1,newdat) )
acc=function(x, y){sum(x==y)/length(x)}
err.glm=acc(p1, glearn)
err.glm

out.glm=classif.glm(glearn~x, data = dat0)
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
  x=mlearn$data,
  y=glearn,
  method="nnet",
  trace=FALSE)
nn
confusionMatrix(nn)
#random forests
rf<-train(
  x=mlearn$data,
  y=glearn,
  method="ranger")
rf
confusionMatrix(rf)

out.glm$max.prob
out.np$max.prob
out.knn$max.prob
out.kernel$max.prob
