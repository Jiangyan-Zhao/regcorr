#' Fit a Regression Model for Pearson Correlation Coefficient
#'
#' Estimate regression coefficients in a covariate-dependent Pearson
#' correlation model for bivariate normal responses using the
#' Newton-Raphson algorithm.
#'
#' @param Y Response matrix with two columns.
#' @param X Covariate matrix including an intercept column.
#' @param betaIni Initial value of regression coefficients.
#' @param link Correlation link function:
#'   \code{"1"} for logistic,
#'   \code{"2"} for tanh.
#'
#' @return A list containing:
#' \describe{
#'   \item{betaCurrent}{Estimated regression coefficients.}
#'   \item{numIter}{Number of Newton-Raphson iterations.}
#'   \item{restart}{Number of restarts.}
#' }
#'
#' @examples
#' set.seed(123)
#'
#' dat <- genDataBN(
#'   numSample = 30,
#'   p = 1,
#'   betaTrue = c(0.5, 0.2),
#'   eta1True = c(0, 0),
#'   eta2True = c(0, 0),
#'   link = "1"
#' )
#'
#' fit <- NRfitBivNormal(
#'   Y = dat$Y,
#'   X = dat$X,
#'   betaIni = c(0, 0),
#'   link = "1"
#' )
#'
#' fit$numIter
#'
#' @export
NRfitBivNormal <- function(Y,X,betaIni,link)
{
  # return new estimate of beta using Newton Raphson method
  #
  # input:
  #        X:   n by p, covariate matrix including first column of ones
  #        Y:   n by 2, paired responses
  #    beta0:   initial estimate of beta


  # initial parameters
  numSample = nrow(Y)
  p = ncol(X)

  TOL=0.01

  residual=lm(Y~X[,-1])$residuals      # residuals of linear reg.
  sigmaHat=sqrt(colMeans(residual^2))  # std of y1 and y2, 1 by 2
  Ytilde=residual%*%diag(1/sigmaHat)   # (y-muHat)/sigma, n by 2
  T1=rowSums(Ytilde^2)                 # n by 1
  T2=apply(Ytilde,1,prod)              # n by 1

  # initial beta
  betaPrev=rep(0,p)   # dummy start
  betaCurrent=betaIni
  numIter=0
  restart=0

  # main loop for Newton-Raphson
  while(sum((betaCurrent-betaPrev)^2)>TOL & numIter<10 & restart<10) # when beta not converge
  {
    numIter=numIter+1
    betaPrev=betaCurrent

    switch(link,
           "1" = { # logistic link
             rho=logistic(X%*%betaCurrent)                       # n by 1, corr. coef.
             a=-rho^2/(1+rho)^2/(1-rho)                          # n by 1
             b=rho*(1+rho^2)/(1+rho)^2/(1-rho)                   # n by 1
             c=rho^2/(1+rho)                                     # n by 1
             u=-rho^2*(2-rho+rho^2)/(1+rho)^3/(1-rho)            # n by 1
             v=-rho*(-1+rho-5*rho^2+rho^3)/(1+rho)^3/(1-rho)     # n by 1
             w=rho^2*(2+rho)*(1-rho)/(1+rho)^2                   # n by 1
           },
           "2" = { # tanh link
             rho=tanh(X%*%betaCurrent)                       # n by 1, corr. coef.
             a=-rho/(1-rho)/(1+rho)                          # n by 1
             b=(1+rho^2)/(1-rho)/(1+rho)                     # n by 1
             c=rho                                           # n by 1
             u=-(1+rho^2)/(1-rho)/(1+rho)                    # n by 1
             v=4*rho/(1-rho)/(1+rho)                         # n by 1
             w=(1-rho)*(1+rho)                               # n by 1
           }
    ) # end of switch

    abc=a*T1+b*T2+c                                     # n by 1
    uvw=u*T1+v*T2+w                                     # n by 1

    #print(uvw)

    score=t(X) %*% abc                                  # p by 1
    H=t(X) %*% diag(as.vector(uvw)) %*% X               # p by p

   # print(H)
   # print(eigen(H))

  if (sqrt(sum(betaCurrent^2))>10 || kappa(H)>10000 ) { # ill posed matrix
    betaPrev=runif(p); betaCurrent=betaIni; numIter=1
    restart=restart+1
  } else {                                            # regular matrix
    invH=solve(H)                                     # p by p, inverse Hessian
    betaCurrent <- betaCurrent-invH%*%score           # update beta
  }  # end of if

} # end of while
#if (restart>1) {print(paste("restart=",restart))}

  return(result=list(betaCurrent=betaCurrent,numIter=numIter,restart=restart))
}
