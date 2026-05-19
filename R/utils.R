#' Logistic function
#'
#' @param x A numeric vector.
#' @return The calculated logistic probability.
#' @export
logistic <- function(x){
  return(1/(1+exp(-x)))
}

#' Generate bivariate binary data
#'
#' @param n Number of rows.
#' @param p 1 by 2 mean vector of bivariate variables.
#' @param rho Correlation of bivariate variables.
#' @return n by 2 matrix of generated binary variables.
#' @importFrom stats rbinom
#' @export
rbinary <- function(n, p, rho) {

  b = rho * sqrt(p[2]*(1-p[2])/p[1]/(1-p[1]))
  Y = matrix(0, n, 2)
  Y[,1] = rbinom(n, 1, p[1])
  Y[,2] = rbinom(n, 1, p[2] + b*(Y[,1]-p[1]))
  return(Y)
}
