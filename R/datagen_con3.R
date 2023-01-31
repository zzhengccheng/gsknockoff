#' This function generate simulated data for gaussian outcomes for K=3
#' @param n1 sample size in population 1
#' @param n2 sample size in population 2
#' @param n3 sample size in population 3
#' @param p number of predictors
#' @param s0 number of true signals
#' @param s1 number of signals shown in population 1 only
#' @param s2 number of signals shown in population 2 only
#' @param s3 number of signals shown in population 3 only
#' @param s12 number of signals shown in population 1 and 2 only
#' @param s13 number of signals shown in population 1 and 3 only
#' @param s23 number of signals shown in population 2 and 3 only
#' @param rho1 correlation parameter for design matrix in population 1
#' @param rho2 correlation parameter for design matrix in population 2
#' @param rho2 correlation parameter for design matrix in population 3
#' @param sigma1 noise parameter for population 1
#' @param sigma2 noise parameter for population 2
#' @param sigma3 noise parameter for population 3
#' @param samesig whether true signals magnitude the same for all data
#' @param effect scale multiplier
#' @return A simulated dataset
#' @export

datagen_con3<-function(n1,n2,n3,p,s0,s1,s2,s3,s12,s13,s23,rho1,rho2,rho3,sigma1,sigma2,sigma3,samesig,scale){

	Sigma1<-Sigma_AR(rho1,p)
	Sigma2<-Sigma_AR(rho2,p)
	Sigma3<-Sigma_AR(rho3,p)
	X1<-mvrnorm(n=n1,mu=rep(0,p),Sigma=Sigma1)
	X2<-mvrnorm(n=n2,mu=rep(0,p),Sigma=Sigma2)
	X3<-mvrnorm(n=n3,mu=rep(0,p),Sigma=Sigma3)
	b0.1<-runif(s0)
	b0.2<-runif(s0)
	b0.3<-runif(s0)
	if (samesig){b0.2<-b0.1;b0.3<-b0.1}
	b1<-runif(s1)
	b2<-runif(s2)
	b3<-runif(s3)
	b12.1<-runif(s12)
	b12.2<-runif(s12)
	b13.1<-runif(s13)
	b13.3<-runif(s13)
	b23.2<-runif(s23)
	b23.3<-runif(s23)
	if (samesig){b12.2<-b12.1;b13.3<-b13.1;b23.3<-b23.2}
	ss0<-sample(c(-1,1),s0,prob=c(0.5,0.5),replace=TRUE)
	ss1<-sample(c(-1,1),s1,prob=c(0.5,0.5),replace=TRUE)
	ss2<-sample(c(-1,1),s2,prob=c(0.5,0.5),replace=TRUE)
	ss3<-sample(c(-1,1),s3,prob=c(0.5,0.5),replace=TRUE)
	ss12<-sample(c(-1,1),s12,prob=c(0.5,0.5),replace=TRUE)
	ss13<-sample(c(-1,1),s13,prob=c(0.5,0.5),replace=TRUE)
	ss23<-sample(c(-1,1),s23,prob=c(0.5,0.5),replace=TRUE)


	beta1<-c(b0.1*ss0,b1*ss1,0*ss2,0*ss3,b12.1*ss12,b13.1*ss13,0*ss23,rep(0,p-s0-s1-s2-s3-s12-s13-s23))*scale
	beta2<-c(b0.2*ss0,0*ss1,b2*ss2,0*ss3,b12.2*ss12,0*ss13,b23.2*ss23,rep(0,p-s0-s1-s2-s3-s12-s13-s23))*scale
	beta3<-c(b0.3*ss0,0*ss1,0*ss2,b3*ss3,0*ss12,b13.3*ss13,b23.3*ss23,rep(0,p-s0-s1-s2-s3-s12-s13-s23))*scale

	y1<-X1%*%beta1+rnorm(n1)*sigma1
	y2<-X2%*%beta2+rnorm(n2)*sigma2
	y3<-X3%*%beta3+rnorm(n3)*sigma3
	data1<-cbind(1,y1,X1)
	data2<-cbind(2,y2,X2)
	data3<-cbind(3,y3,X3)
	data<-data.frame(rbind(data1,data2,data3))
	names(data)=c("C","Y",paste("X",1:p,sep=""))
	data
}