#' Generate Correlated Binary Data
#'
#' Generate bivariate binary random variables with specified marginal
#' probabilities and correlation.
#'
#' The implementation follows Qaqish (2003).
#'
#' @param n Number of observations.
#' @param p A length-2 vector of marginal probabilities.
#' @param rho Correlation coefficient.
#'
#' @return An \code{n x 2} matrix of binary random variables.
#'
#' @references
#' Qaqish, B. F. (2003). A family of multivariate binary distributions for
#' simulating correlated binary variables with specified marginal means and
#' correlations. \emph{Biometrika}, 90(2), 455--463.
#'
#' @examples
#' set.seed(123)
#' y <- rbinary(
#'   n = 10,
#'   p = c(0.4, 0.5),
#'   rho = 0.2
#' )
#'
#' dim(y)
#'
#' @export
rbinary <-function(n,p,rho)
{
b=rho*sqrt(p[2]*(1-p[2])/p[1]/(1-p[1]))
Y=matrix(0,n,2)
Y[,1]=rbinom(n,1,p[1])
Y[,2]=rbinom(n,1,p[2]+b*(Y[,1]-p[1]))
return(Y)
}
