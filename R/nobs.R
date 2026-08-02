#' Number of Observations
#'
#' Returns the number of observations used to fit a \code{regcorr} model.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The number of observations, a single integer.
#'
#' @seealso \code{\link{regcorr}} for model fitting.
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
#' nobs(fit)
#'
#' @name nobs
#' @importFrom stats nobs
#' @export
nobs.regcorr <- function(object, ...) {
  return(nrow(object$y))
}
