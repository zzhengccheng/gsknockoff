#' This function implement our proposed Knockoff method for testing union null 
#' @param Xlist list of length K for predictor matrices from each experiment
#' @param Ylist list of length K for outcomes from each experiment
#' @param familylist vector of length K for type of outcomes from each experiment (Gaussian, Binomial or Cox)
#' @param Kmethod vecttor of length K for knockoff construction methods for each experiment (fixed or second order)
#' @param Zmethod vecttor of length K for Z statistics calculation methods for each experiment. Current implemented method is glmnet   
#' @param Zfun User defined function on how to construct Z statistics from a list of length K for Z and a list of length K for Ztilde (required if select Zmethod="Other")
#' @param q FDR level for control
#' @param offset 0 for Knockoff and 1 for Knockoff+
#' @param Wmethod Current implemented method for general K is Diff, method="Other" can be used for other user defined function 
#' @param Wfun User defined function on how to construct W statistics from a list of length K for Z and a list of length K for Ztilde (required if select Wmethod="Other")
#' @return A list of selected variables
#' @export
myest_sim<-function(Xlist,Ylist,familylist,Kmethod,Zmethod,Wmethod,q,offset,Zfun=NULL,Wfun=NULL){
	X<-list()
	Z<-list()
	Xc<-list()
	XKnock<-list()
	Ztilde<-list()
	p=ncol(Xlist[[1]])
	K=length(Ylist)
	for (k in 1:K){
		X[[k]]=as.matrix(Xlist[[k]])
		Xc[[k]]<-X[[k]]%*%diag(1/sqrt(colSums(X[[k]]^2)))		
		if (Kmethod[k]=="fixed"){
			XKnock[[k]]<-create.fixed(Xc[[k]])$Xk
		}
		if (Kmethod[k]=="second"){
			XKnock[[k]]<-create.second_order(Xc[[k]],shrink=T)
		}
	}
	for (k in 1:K){
		y<-Ylist[[k]]
		if (Zmethod[k]=="glmnet"){
			if (familylist[k]=="binomial"){
				fit=cv.glmnet(as.matrix(cbind(Xc[[k]],XKnock[[k]])),c(y),family=familylist[k],standardize=FALSE,intercept=TRUE)
			}else{
				fit=cv.glmnet(as.matrix(cbind(Xc[[k]],XKnock[[k]])),c(y),family=familylist[k],standardize=FALSE,intercept=FALSE)
			}
			if (familylist[k]=="cox"){
				Z[[k]]=abs(coef(fit,fit$lambda.min)[(1:p)])	
				Ztilde[[k]]=abs(coef(fit,fit$lambda.min)[p+(1:p)])
			}else{
				Z[[k]]=abs(coef(fit)[1+(1:p)])	
				Ztilde[[k]]=abs(coef(fit)[1+p+(1:p)])
			}
		}	
		if (Zmethod[k]=="Other"){
			tmp=Zfun(y, Xc[[k]], XKnock[[k]], familylist[k])
			Z[[k]]=tmp[1:p]
			Ztilde[[k]]=tmp[p+(1:p)]
		}		
	}
	if (Wmethod=="Diff"){
		W=rep(1,p)
		for (k in 1:K){
			W=W*(Z[[k]]-Ztilde[[k]])
		}
	}
	if (Wmethod=="Other"){
		W=Wfun(Z,Ztilde)
	}
	mythred=knockoff.threshold(W,fdr=q,offset=offset)
	myselect=which(W>=mythred)
	return(myselect)
}