#' Logistic function
#'
#' @param x A numeric vector.
#' @return The calculated logistic probability.
logistic <- function(x){
  return(stats::plogis(x))
}

#' Generate bivariate binary data
#'
#' @param n Number of rows.
#' @param p 1 by 2 mean vector of bivariate variables.
#' @param rho Correlation of bivariate variables.
#' @return n by 2 matrix of generated binary variables.
#' @importFrom stats rbinom
rbinary <- function(n, p, rho) {

  b = rho * sqrt(p[2]*(1-p[2])/p[1]/(1-p[1]))
  Y = matrix(0, n, 2)
  Y[,1] = rbinom(n, 1, p[1])
  # conditional success probability; clamp to [0,1] against floating-point
  # violations of the feasible range
  p2cond = pmin(pmax(p[2] + b*(Y[,1]-p[1]), 0), 1)
  Y[,2] = rbinom(n, 1, p2cond)
  return(Y)
}
