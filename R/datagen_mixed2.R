#' This function generate simulated data for mixed outcomes for K=2 (one binary, one gaussian)
#' @param n1 sample size in population 1
#' @param n2 sample size in population 2
#' @param p number of predictors
#' @param s0 number of true signals
#' @param s1 number of signals shown in population 1 only
#' @param s2 number of signals shown in population 2 only
#' @param rho1 correlation parameter for design matrix in population 1
#' @param rho2 correlation parameter for design matrix in population 2
#' @param sigma1 noise parameter for population 1
#' @param sigma2 noise parameter for population 2
#' @param samesig whether true signals magnitude the same for all data
#' @param effect scale multiplier
#' @return A simulated dataset
#' @export
datagen_mixed2<-function(n1,n2,p,s0,s1,s2,rho1,rho2,sigma1,sigma2,samesig,scale){
	Sigma1<-Sigma_AR(rho1,p)
	Sigma2<-Sigma_AR(rho2,p)
	X1<-mvrnorm(n=n1,mu=rep(0,p),Sigma=Sigma1)
	X2<-mvrnorm(n=n2,mu=rep(0,p),Sigma=Sigma2)
	b0.1<-runif(s0)
      b0.2<-runif(s0)
      if (samesig){b0.2<-b0.1}
	b1<-runif(s1)
	b2<-runif(s2)
	ss0<-sample(c(-1,1),s0,prob=c(0.5,0.5),replace=TRUE)
	ss1<-sample(c(-1,1),s1,prob=c(0.5,0.5),replace=TRUE)
	ss2<-sample(c(-1,1),s2,prob=c(0.5,0.5),replace=TRUE)
	beta1<-c(b0.1*ss0,b1*ss1,rep(0,p-s0-s1))*scale
 	beta2<-c(b0.2*ss0,rep(0,s1),b2*ss2,rep(0,p-s0-s1-s2))*scale
	y1<-X1%*%beta1+rnorm(n1)*sigma1
	y2<-X2%*%beta2+rnorm(n2)*sigma2
	####dichotomize y2 (same as probit model)
	y2<-as.numeric(y2>0)
	data1<-cbind(1,y1,X1)
	data2<-cbind(2,y2,X2)
	data<-data.frame(rbind(data1,data2))
	names(data)=c("C","Y",paste("X",1:p,sep=""))
	data
}