#' @examples
#' set.seed(123)
#' dat <- genDataBB(
#'   numSample = 20,
#'   p = 2,
#'   betaTrue = c(0.2, 0, 0),
#'   eta1True = c(0, 0, 0),
#'   eta2True = c(0, 0, 0),
#'   link = "1"
#' )
#' names(dat)
genDataBB<-function(numSample,p,betaTrue,eta1True,eta2True,link)
{ # return data generated from biv Bernoulli

  # generate X
  X=cbind(rep(1,numSample),matrix(runif(numSample*p),numSample))

  # generate rho
  switch(link, # 1: logistic; 2:tanh
         "1" = {rho=as.matrix(logistic(X%*%betaTrue),ncol=1) },
         "2" = {rho=as.matrix(tanh(X%*%betaTrue),ncol=1) }
  )        # end of switch

  # generate true prob of success and true corr. coef.
  p1=as.matrix(logistic(X%*%eta1True),ncol=1)   # marginal success prob for y1, numSample by 1
  p2=as.matrix(logistic(X%*%eta2True),ncol=1)   # marginal success prob for y2, numSample by 1

  # validity check for the legal relationship between these p1, p2 and rho
  psi1=sqrt(p1/(1-p1)); psi2=sqrt(p2/(1-p2))
  validID=apply(cbind(-psi1*psi2,-1/(psi1*psi2)),2,max)<=rho
  validID=validID & rho<=apply(cbind(psi1/psi2,psi2/psi1),2,min)       # numSample by 1 of logic values

  # filter out some invalid p1, p2, rho
  p1=p1[validID]; p2=p2[validID]
  rho=rho[validID]                      # in this example, they are all valid.
  X=X[validID,]
  numSample=length(p1)                  # re-count numSample

  # generate corresponding (Y1,Y2)
  Y=matrix(0,numSample,2)               # define initial Y, numSample by 2
  for (iSample in 1:numSample)
  {Y[iSample,]=rbinary(1,c(p1[iSample],p2[iSample]),rho[iSample]) }       # each pair of (y1, y2)

  return(list(X=X,Y=Y,rho=rho))
}
