#' Extract the Variance-Covariance Matrix
#'
#' Extracts the variance-covariance matrix of the estimated coefficients
#' from a fitted \code{regcorr} model. The matrix is computed from
#' bootstrap replications of the coefficients; if the model was fitted with
#' \code{nboot = 0}, the matrix contains \code{NA} values.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A square matrix with one row and column per coefficient.
#'
#' @seealso \code{\link{regcorr}} for model fitting, especially the
#'   \code{nboot} argument; \code{\link{summary.regcorr}} for standard
#'   errors and p-values; \code{\link{coef.regcorr}} for the coefficients.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' x <- runif(n)
#' rho <- plogis(0.5 + 0.2 * x)
#' z1 <- rnorm(n)
#' y1 <- z1
#' y2 <- z1 * rho + rnorm(n) * sqrt(1 - rho^2)
#' fit <- regcorr(cbind(y1, y2) ~ x, data = data.frame(y1, y2, x), nboot = 20)
#' vcov(fit)
#'
#' @name vcov
#' @importFrom stats vcov
#' @export
vcov.regcorr <- function(object, ...) {
  return(object$vcov)
}
