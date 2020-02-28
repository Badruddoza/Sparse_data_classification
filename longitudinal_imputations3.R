#install.packages("longitudinalData")
rm(list=ls())
set.seed(123)
setwd("E:/2. Other research/1. Classification of Gaussian Process/actual data")
#install.packages("dplyr")
require(dplyr);require(longitudinalData);require(tidyverse);require(tidyr);
require(ggplot2);require(fda.usc);library(knitr);library(kableExtra);library(caret)
#install.packages("magick");#install.packages("webshot");#webshot::install_phantomjs()
#install.packages('processx')

################################################################# load data set
dat=readRDS("Adrenal-EMRE-dynamic.RDS")
#separate out the working data and rename variables
dat0<-cbind.data.frame(mal=as.integer(dat$malignancy),
                       id=as.factor(dat$lesionID),
                       time=as.integer(dat$time),
                       mean=as.numeric(dat$calc_mean))

################################################################# summary statistics
dat_sum<-rbind.data.frame(cbind("Malignancy","Dummy variable",round(mean(dat0$mal),3),round(sd(dat0$mal)),round(min(dat0$mal),3),round(max(dat0$mal),3)),
  cbind("Mean index","Numeric",round(mean(dat0$mean),3),round(sd(dat0$mean),3),round(min(dat0$mean),3),round(max(dat0$mean),3)),
  cbind("Time of scan","Integer",round(mean(dat0$time),3),round(sd(dat0$time),3),min(dat0$time),max(dat0$time)),
  cbind("Patient ID","Character","Count=",length(unique(dat0$id)),"Sample=",dim(dat0)[1]))
colnames(dat_sum)<-c("Variable","Type","Mean","SD","Min","Max")

kable(dat_sum, caption="Table 1: Summary Statistics") %>%
  kable_styling(bootstrap_options=c("striped","hover")) %>%
  save_kable("sumtable.png")

################################################################# plot the data set
ggplot(data=dat0,aes(x=time,y=mean,color=as.factor(mal)))+
  geom_point()+
  geom_line(aes(group=id), alpha = 0.2)+
  ggtitle("Figure 1: Scan results over time")+
  xlab("Time of scan")+
  ylab("Mean Index")+
  labs(color="Malignancy")+
  theme_bw()
ggsave("imp.png",width=4,height=4)

################################################################# split into training and testing sample
dat1<-spread(dat0,time,mean)
dat1<-dat1 %>% dplyr::select(id,names(dat1)[which(names(dat1)!="id")]) # make id the first variable
trainIndex<-sample(seq_len(dim(dat1)[1]),size=floor(0.9*nrow(dat1)))

########################################################## THE LONG LOOP ##################################
if(file.exists("results.csv")){
  #file.remove("results.csv")
  print("peo")
}
i=1
j="copyMean.bisector"
method_list<-c("linearInterpol.locf", "linearInterpol.local", "linearInterpol.bisector", 
               "trajMean", "trajMedian", "trajHotDeck", "locf", "nocb", 
               "copyMean.locf", "copyMean.local", "copyMean.bisector")
#simplest_method_list<-c("crossMean", "crossMedian", "crossHotDeck")
#not_running<-c("linearInterpol.global","copyMean.global")
# locf=last observation carried forward, nocb=next observation carried before

for(j in method_list){
  
print(paste0("Starting method ",j," at ",Sys.time()))
keep<-c("dat0","dat1","trainIndex","i","j")
rm(list=ls()[!ls() %in% keep])

dat1_train<-dat1[trainIndex,];dat1_test<-dat1[-trainIndex,]
dat1_train_m<-dat1_train[dat1_train$mal==1,];dat1_train_b<-dat1_train[dat1_train$mal==0,]
dat1_test_n<-dat1_test %>% mutate(mal=NA) # make the response of test data missing
dat1_m<-rbind(dat1_train_m,dat1_test_n);dat1_b<-rbind(dat1_train_b,dat1_test_n)

################################################################# imputation
dat1_m[,2:dim(dat1_m)[2]]<-imputation(as.matrix(dat1_m[,2:dim(dat1_m)[2]]),j)
dat1_b[,2:dim(dat1_b)[2]]<-imputation(as.matrix(dat1_b[,2:dim(dat1_b)[2]]),j)

################################################################# imputation performance by malignant and benign sample in test data
# first identify which part of the sample was test
dat1_m_i<-dat1_m[which(dat1_m$id %in% dat1_test$id),];dat1_b_i<-dat1_b[which(dat1_b$id %in% dat1_test$id),]
# make the response missing because it got imputed meanwhile
dat1_m_i<-dat1_m_i %>% mutate(mal=NA);dat1_b_i<-dat1_b_i %>% mutate(mal=NA)
# count nonmissing columns in each low
dat1_test1<-dat1_test %>% mutate(rs=dim(dat1_test)[2]-2-rowSums(is.na(dat1_test[,3:dim(dat1_test)[2]])))
adat<-bind_rows(cbind(dat1_m_i,m="Malignant" ),
                cbind(dat1_b_i,m="Benign"),
                cbind(dat1_test,m="Actual"))
# make the data long for graphs
adat1<-gather(adat,"time","mean",3:(dim(adat)[2]-1));
adat1<-adat1 %>% arrange(m,id,time)
adat2<-spread(adat1,"m","mean")
adat2<-adat2 %>% arrange(id,time) %>%
  group_by(id,time) %>%
  summarize_at(c("mal","Actual","Benign","Malignant"),mean,na.rm=TRUE)
adat2$Actual[is.nan(adat2$Actual)==TRUE]<-NA  #replace Nan with NA

########### Plotting selected points to see the performance of imputation
# keep only 3 lesions with 1,2, and 3 observations
shortlist<-dat1_test1 %>% filter(rs<=3) %>%  dplyr::select(id,rs) %>%  
  arrange(rs) %>%  mutate(rs1=ifelse(rs!=lag(rs),1,0)) %>%
  mutate(rs1=ifelse(is.na(rs1)==TRUE,1,rs1)) %>%  filter(rs1==1)
gdat<-adat2 %>% filter(id %in% shortlist$id) %>%
  rename(Benign_based=Benign) %>%
  rename(Malignant_based=Malignant)
gdat<-gather(gdat,"Imputation","mean",4:(dim(gdat)[2]))
ggplot()+
  geom_line(data=gdat[gdat$Imputation!="Actual",], aes( as.numeric(time),as.numeric(mean),group=id,color=Imputation),lwd=1)+
  geom_point(data=gdat[gdat$Imputation=="Actual",],aes(as.numeric(time),as.numeric(mean),shape=as.factor(mal)),size=3)+
  labs(title=paste0("Figure ",i+1,": Imputed by ",j),x="Time of scan", y="Mean index")+theme_bw()+
  theme(plot.title=element_text(size=12))+
  scale_shape(name="Malignancy")

ggsave(paste0("p_",i,".png"),width=4,height=4)
print("graph saved")



############################# choosing the best method
if(sum((adat2$Actual-adat2$Benign)^2,na.rm=TRUE)>=sum((adat2$Actual-adat2$Malignant)^2,na.rm=TRUE)){
  print("Malignant performs better")
  data<-merge(dat1_m,dat1_b[which(dat1_b$id %in% dat1_train_b$id)],all=TRUE)
  }else{
    print("Benign performs better")
    data<-merge(dat1_b,dat1_m[which(dat1_m$id %in% dat1_train_m$id)],all=TRUE)
    }
data<-data %>% mutate(test=ifelse(mal!=0 & mal!=1,1,0)) %>% 
  dplyr::select(test,names(data)[which(names(data)!="test")]) %>% arrange(test,id)
data[,3:dim(data)[2]]<-imputation(as.matrix(data[,3:dim(data)[2]]),j)
data$mal[data$test==1]<-NA
data$mal<-as.factor(data$mal)

######################################################### Analysis
data_train<-data %>% filter(test==0) 
data_test<-data %>% filter(test==1)

#using LongitudinalData package
mdata1<-data_train[,4:dim(data_train)[2]]
mlearn1=fdata(mdata1)
glearn=data_train$mal
dataf<-data.frame(glearn) 
data_train_l=list("df"=dataf,"x"=mlearn1)

out.glm=classif.glm(glearn~x, data = data_train_l)
out.np=classif.np(glearn,mlearn1)
out.knn=classif.knn(glearn,mlearn1,knn=c(3,5,7))
out.kernel=classif.kernel(glearn,mlearn1)

mdata2<-data_test[,4:dim(data_test)[2]]
mlearn2=fdata(mdata2)
data_test_l<-list("x"=mlearn2)

glm1<-predict.classif(out.glm,data_test_l)
np1<-predict.classif(out.np,mlearn2)
knn1<-predict.classif(out.knn,mlearn2)
kernel1<-predict.classif(out.kernel,mlearn2)
print("LongitudinalData completed")

#using machine learning
#neural network
nn<-train(x=mdata1,y=data_train$mal,method="nnet",trace=FALSE)
nn1<-predict(nn,mdata2)

#random forests
rf<-train(x=mdata1,y=data_train$mal,method="ranger")
rf1<-predict(rf,mdata2)

#extreme gradient boosting
gb<-train(x=mdata1,y=data_train$mal,method="xgbTree")
gb1<-predict(gb,mdata2)
print("Machine Learning completed")


#collect predictions
predicted<-merge(dat1_test,cbind.data.frame(id=data_test$id,
                                            GLM=glm1,
                                            NP=np1,
                                            KNN=knn1,
                                            Kernel=kernel1,
                                            NNet=nn1,
                                            RF=rf1,
                                            XGB=gb1),all=TRUE)
methods<-colnames(predicted)[(names(predicted) %in% names(data_test))==FALSE]
predicted1<-predicted %>% dplyr::select(mal,methods)
predicted1<-predicted1 %>% mutate_all(as.character)

# overall performance
reslist0<-c()
for(v in names(predicted1)[-which(names(predicted1)=="mal")]){
  predicted2<-predicted1 %>% dplyr::select(mal,v)
  tt<-table(predicted2)
  a=0;d=0;c=0;b=0
  tryCatch({a<-tt[2,2]}, error=function(e){});
  tryCatch({b<-tt[1,2]}, error=function(e){});
  tryCatch({c<-tt[2,1]}, error=function(e){});
  tryCatch({d<-tt[1,1]}, error=function(e){});
  tn<-d/(sum(tt));tp<-a/(sum(tt));fp<-b/sum(tt);fn<-c/sum(tt)
  tpr<-tp/(tp+fn);tnr=tn/(tn+fp);acc<-(tp+tn)/(tp+tn+fp+fn);
  nn=(a+b+c+d);po=(a+d)/nn;py=((a+b)/nn)*((a+c)/nn);pn=((c+d)/nn)*((b+d)/nn)
  pe=py+pn;kappa=(po-pe)/(1-pe);
  reslist0<-c(reslist0,v,round(tpr,4),round(tnr,4),round(acc,4),round(kappa,4) )
}
reslist1<-data.frame(matrix(c(reslist0),nrow=length(methods),byrow=TRUE)) %>%
  arrange(-as.numeric(X5)) %>% mutate(Scan="Overall") %>%
  rename(Model=X1) %>% rename(TPR=X2) %>% rename(TNR=X3) %>% rename(Accuracy=X4) %>% rename(Kappa=X5)
#reslist1

# performance with one multiple scans only
tvars<-colnames(predicted)[(names(predicted) %in% names(data_test))==TRUE]
tvars<-tvars[-c(1,2)]
predicted3<-predicted %>%
  mutate(Scan=dim(predicted[,tvars])[2]-rowSums(is.na(predicted[,tvars]))) %>%
  dplyr::select(mal,methods,Scan) %>%
  filter(Scan>1) %>% dplyr::select(-Scan)
predicted3<-predicted3 %>% mutate_all(as.character)

reslist0<-c()
for(v in names(predicted3)[-which(names(predicted3)=="mal")]){
  predicted4<-predicted3 %>% dplyr::select(mal,v)
  tt<-table(predicted4)
  a=0;d=0;c=0;b=0
  tryCatch({a<-tt[2,2]}, error=function(e){});
  tryCatch({b<-tt[1,2]}, error=function(e){});
  tryCatch({c<-tt[2,1]}, error=function(e){});
  tryCatch({d<-tt[1,1]}, error=function(e){});
  tn<-d/(sum(tt));tp<-a/(sum(tt));fp<-b/sum(tt);fn<-c/sum(tt)
  tpr<-tp/(tp+fn);tnr=tn/(tn+fp);acc<-(tp+tn)/(tp+tn+fp+fn)
  nn=(a+b+c+d);po=(a+d)/nn;py=((a+b)/nn)*((a+c)/nn);pn=((c+d)/nn)*((b+d)/nn)
  pe=py+pn;kappa=(po-pe)/(1-pe);
  reslist0<-c(reslist0,v,round(tpr,4),round(tnr,4),round(acc,4),round(kappa,4) )
}
reslist2<-data.frame(matrix(c(reslist0),nrow=length(methods),byrow=TRUE)) %>%
  arrange(-as.numeric(X5)) %>% mutate(Scan="Multiple") %>%
  rename(Model=X1) %>% rename(TPR=X2) %>% rename(TNR=X3) %>% rename(Accuracy=X4) %>% rename(Kappa=X5)
#reslist2

# Combine all results together
reslist<-cbind(Impuation=j,rbind.data.frame(reslist1,reslist2))

kable(reslist, caption=paste0("Table ",i+1,": Model Performance") ) %>%
  kable_styling(bootstrap_options=c("striped","hover")) %>%
  save_kable(paste0("performance_",i,".png"))
write.table(reslist,"results.csv",append=T,sep=",",row.names = F,col.names = F)


i=i+1
print(paste0("Done with method ",j," at ",Sys.time() ) )
}