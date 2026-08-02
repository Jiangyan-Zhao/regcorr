#' Generate bivariate normal data
#'
#' Simulate bivariate normal responses with a covariate-dependent
#' Pearson correlation coefficient.
#'
#' @param numSample Number of observations.
#' @param p Number of covariates.
#' @param betaTrue True regression coefficients for the correlation model.
#' @param eta1True Included for interface consistency.
#' @param eta2True Included for interface consistency.
#' @param link Correlation link function: "1" for logistic, "2" for tanh.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{Covariate matrix.}
#'   \item{Y}{Response matrix with two columns.}
#'   \item{rho}{True correlation coefficients.}
#' }
#'
#' @examples
#' set.seed(123)
#' dat <- genDataBN(
#'   numSample = 20,
#'   p = 2,
#'   betaTrue = c(0.2, 0, 0),
#'   eta1True = c(0, 0, 0),
#'   eta2True = c(0, 0, 0),
#'   link = "1"
#' )
#' names(dat)
#'
#' dim(dat$X)
#' dim(dat$Y)
#'
#' @export
genDataBN<-function(numSample,p,betaTrue,eta1True,eta2True,link)
{ # return data generated from biv normal

  # generate X, numSample by numCova+1, all uniform
  X=cbind(rep(1,numSample),matrix(runif(numSample*p),numSample))

    switch(link, # 1: logistic; 2:tanh
         "1" = {rho=as.matrix(logistic(X%*%betaTrue),ncol=1) },
         "2" = {rho=as.matrix(tanh(X%*%betaTrue),ncol=1) }
  ) # end of switch

  sigma1=1; sigma2=1
  # generate Y with eta1 and eta2 both zero vector
  Y=matrix(0,numSample,2)
  z1=rnorm(numSample)
  Y[,1]=sigma1*z1                                       # Y[,1]~N(0,sigma1^2)
  Y[,2]=sigma2*(z1*rho+rnorm(numSample)*sqrt(1-rho^2))  # Y[,2]~N(0,sigma2^2)

  return(list(X=X,Y=Y,rho=rho))
}
