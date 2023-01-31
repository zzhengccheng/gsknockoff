#' This function implements the Pooling Knockoff method
#' @param Xlist list of length K for predictor matrices from each experiment
#' @param Ylist list of length K for outcomes from each experiment
#' @param familylist list of length K for type of outcomes from each experiment (Gaussian, Binomial or Cox)
#' @param Kmethod list of length K for knockoff construction methods for each experiment (fixed or second order)
#' @param Zmethod Current implemented method is glmnet  
#' @param Zfun User defined function on how to construct Z statistics from a list of length K for Z and a list of length K for Ztilde (required if select Zmethod="Other")
#' @param q FDR level for control
#' @param offset 0 for Knockoff and 1 for Knockoff+
#' @return A list of selected variables
#' @export

myest_stack<-function(Xlist,Ylist,familylist,Kmethod,Zmethod,q,offset,Zfun=NULL,IDorder=NULL){
	####check whether familylist, Kmethod, Zmethod are the same
	samefamily=as.numeric(length(unique(familylist))==1)
	sameKmethod=as.numeric(length(unique(Kmethod))==1)
	sameZmethod=as.numeric(length(unique(Zmethod))==1)
	
	if (samefamily & sameKmethod & sameZmethod){
		p=ncol(Xlist[[1]])
		K=length(Ylist)

		X<-as.matrix(Xlist[[1]])
		YY<-Ylist[[1]]
		if (K>1){
			for (k in 2:K){
				X=rbind(X,as.matrix(Xlist[[k]]))
				YY=c(YY,Ylist[[k]])
			}
		}

		if (!is.null(IDorder)){		
			X[IDorder,]<-X
			if (familylist[1]=="cox"){
				YY[IDorder,]<-YY
			}else{
				YY[IDorder]<-YY
			}
		}
		Xc<-X%*%diag(1/sqrt(colSums(X^2)))
		if (Kmethod[1]=="fixed"){
			XKnock<-create.fixed(Xc)$Xk
		}
		if (Kmethod[1]=="second"){
			XKnock<-create.second_order(Xc,shrink=T)
		}
		if (Zmethod[1]=="glmnet"){
			if (familylist[1]=="binomial"){
				fit=cv.glmnet(as.matrix(cbind(Xc,XKnock)),c(YY),family=familylist[1],standardize=FALSE,intercept=TRUE)
			}else{
				fit=cv.glmnet(as.matrix(cbind(Xc,XKnock)),c(YY),family=familylist[1],standardize=FALSE,intercept=FALSE)
			}
			if (familylist[1]=="cox"){
				Z=abs(coef(fit,fit$lambda.min)[(1:p)])	
				Ztilde=abs(coef(fit,fit$lambda.min)[p+(1:p)])
			}else{
				Z=abs(coef(fit)[1+(1:p)])	
				Ztilde=abs(coef(fit)[1+p+(1:p)])
			}
		}	
		if (Zmethod[1]=="Other"){
			tmp=Zfun(y, Xc, XKnock, familylist[1])
			Z=tmp[1:p]
			Ztilde=tmp[p+(1:p)]
		}				
		W=Z-Ztilde
		mythred=knockoff.threshold(W,fdr=q,offset=offset)
		myselect=which(W>=mythred)
		return(myselect)
	}
}