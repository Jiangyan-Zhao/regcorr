#' Extract Model Coefficients
#'
#' Extracts the estimated correlation-link coefficients of a fitted
#' \code{regcorr} model. The coefficients are on the scale of the link
#' function (logistic or tanh); apply the corresponding inverse link
#' (e.g., \code{plogis}) to interpret them as correlations.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named numeric vector of coefficients.
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{vcov.regcorr}} for the variance-covariance matrix;
#'   \code{\link{predict.regcorr}} for predicted correlations.
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
#' coef(fit)
#'
#' @name coef
#' @importFrom stats coef
#' @export
coef.regcorr <- function(object, ...) {
  return(object$coefficients)
}
