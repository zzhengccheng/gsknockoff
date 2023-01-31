####example code
library(devtools)
devtools::install_github("zzhengccheng/gsknockoff")
library(gsknockoff)
library(MASS)
library(knockoff)
library(glmnet)

####K=2 continuous
data=datagen_con2(n1=1000,n2=1000,p=200,s0=40,s1=40,s2=40,rho1=0.4,rho2=0.6,sigma1=1,sigma2=2,samesig=0,scale=2)
Xlist=list(data[which(data$C==1),-c(1:2)],data[which(data$C==2),-c(1:2)])
Ylist=list(data$Y[which(data$C==1)],data$Y[which(data$C==2)])
familylist=c("gaussian","gaussian")
Kmethod=c("fixed","fixed")
Zmethod=c("glmnet","glmnet")
Wmethod="Diff"
myselect=myest_sim(Xlist,Ylist,familylist,Kmethod,Zmethod,Wmethod,q=0.2,offset=1)
myselect_stack=myest_stack(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
myselect_int=myest_int(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
list(myselect,myselect_stack,myselect_int)

####K=2 binary
data=datagen_binary2(n1=1000,n2=1000,p=200,s0=40,s1=40,s2=40,rho1=0.4,rho2=0.6,sigma1=1,sigma2=-1,samesig=0,scale=2)
Xlist=list(data[which(data$C==1),-c(1:2)],data[which(data$C==2),-c(1:2)])
Ylist=list(data$Y[which(data$C==1)],data$Y[which(data$C==2)])
familylist=c("binomial","binomial")
Kmethod=c("second","second")
Zmethod=c("glmnet","glmnet")
Wmethod="Diff"
myselect=myest_sim(Xlist,Ylist,familylist,Kmethod,Zmethod,Wmethod,q=0.2,offset=1)
myselect_stack=myest_stack(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
myselect_int=myest_int(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
list(myselect,myselect_stack,myselect_int)


####K=2 mixed
data=datagen_mixed2(n1=1000,n2=1000,p=200,s0=40,s1=40,s2=40,rho1=0.4,rho2=0.6,sigma1=1,sigma2=2,samesig=0,scale=2)
Xlist=list(data[which(data$C==1),-c(1:2)],data[which(data$C==2),-c(1:2)])
Ylist=list(data$Y[which(data$C==1)],data$Y[which(data$C==2)])
familylist=c("gaussian","binomial")
Kmethod=c("fixed","second")
Zmethod=c("glmnet","glmnet")
Wmethod="Diff"
myselect=myest_sim(Xlist,Ylist,familylist,Kmethod,Zmethod,Wmethod,q=0.2,offset=1)
myselect_stack=myest_stack(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
myselect_int=myest_int(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
list(myselect,myselect_stack,myselect_int)



####K=3 continuous
data=datagen_con3(n1=1000,n2=1000,n3=1000,p=200,s0=40,s1=10,s2=10,s3=10,s12=10,s13=10,s23=10,rho1=0.4,rho2=0.5,rho3=0.6,sigma1=1,sigma2=2,sigma3=1.5,samesig=0,scale=2)
Xlist=list(data[which(data$C==1),-c(1:2)],data[which(data$C==2),-c(1:2)])
Ylist=list(data$Y[which(data$C==1)],data$Y[which(data$C==2)])
familylist=c("gaussian","gaussian","gaussian")
Kmethod=c("fixed","fixed","fixed")
Zmethod=c("glmnet","glmnet","glmnet")
Wmethod="Diff"
myselect=myest_sim(Xlist,Ylist,familylist,Kmethod,Zmethod,Wmethod,q=0.2,offset=1)
myselect_stack=myest_stack(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
myselect_int=myest_int(Xlist,Ylist,familylist,Kmethod,Zmethod,q=0.2,offset=1)
list(myselect,myselect_stack,myselect_int)


