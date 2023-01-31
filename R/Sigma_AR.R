#' This function generates the AR variance covariance matrix with diagonal 1 and correlation rho
#' @param rho correlation for adjacent variables
#' @param p dimension of the matrix
#' @return A variance-covariance matrix
#' @export
Sigma_AR<-function(rho,p){
	rho^abs(outer(1:p,1:p,"-"))
}
