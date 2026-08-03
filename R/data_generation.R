#' Generate data from bivariate normal
#'
#' @param numSample Sample size.
#' @param p Number of covariates.
#' @param betaTrue True beta vector.
#' @param eta1True True eta1 vector.
#' @param eta2True True eta2 vector.
#' @param link Link function indicator ("1" = logistic; "2" = tanh).
#' @return A list containing X, Y, and rho.
#' @importFrom stats runif rnorm
genDataBN <- function(numSample, p, betaTrue, eta1True, eta2True, link) {
  X = cbind(rep(1,numSample), matrix(runif(numSample*p), numSample))
  switch(link,
         "1" = {rho = as.matrix(logistic(X %*% betaTrue), ncol=1)},
         "2" = {rho = as.matrix(tanh(X %*% betaTrue), ncol=1)}
  )
  sigma1=1; sigma2=1
  Y = matrix(0, numSample, 2)
  z1 = rnorm(numSample)
  Y[,1] = sigma1 * z1
  Y[,2] = sigma2 * (z1*rho + rnorm(numSample)*sqrt(1-rho^2))
  return(list(X=X, Y=Y, rho=rho))
}

#' Generate data from bivariate Bernoulli
#'
#' @param numSample Sample size.
#' @param p Number of covariates.
#' @param betaTrue True beta vector.
#' @param eta1True True eta1 vector.
#' @param eta2True True eta2 vector.
#' @param link Link function indicator ("1" = logistic; "2" = tanh).
#' @return A list containing X, Y, and rho.
#' @importFrom stats runif
genDataBB <- function(numSample, p, betaTrue, eta1True, eta2True, link) {
  X = cbind(rep(1,numSample), matrix(runif(numSample*p), numSample))
  switch(link,
         "1" = {rho = as.matrix(logistic(X %*% betaTrue), ncol=1)},
         "2" = {rho = as.matrix(tanh(X %*% betaTrue), ncol=1)}
  )
  p1 = as.matrix(logistic(X %*% eta1True), ncol=1)
  p2 = as.matrix(logistic(X %*% eta2True), ncol=1)

  # feasible range of rho given the marginal probabilities, applied row-wise
  # (Dufera et al., 2023, Section 6; Qaqish, 2003):
  #   max(-psi1*psi2, -(psi1*psi2)^-1) <= rho <= min(psi1/psi2, psi2/psi1)
  psi1 = sqrt(p1/(1-p1)); psi2 = sqrt(p2/(1-p2))
  validID = pmax(-psi1*psi2, -1/(psi1*psi2)) <= rho &
    rho <= pmin(psi1/psi2, psi2/psi1)

  p1 = p1[validID]; p2 = p2[validID]; rho = rho[validID]
  X = X[validID,]
  numSample = length(p1)

  Y = matrix(0, numSample, 2)
  for (iSample in 1:numSample) {
    Y[iSample,] = rbinary(1, c(p1[iSample], p2[iSample]), rho[iSample])
  }
  return(list(X=X, Y=Y, rho=rho))
}
