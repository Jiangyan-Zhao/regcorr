#' Extract Fitted Correlations
#'
#' Extracts the fitted Pearson correlation coefficients for the training
#' data from a fitted \code{regcorr} model.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of fitted correlations.
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{predict.regcorr}} for predicting at new covariate values;
#'   \code{\link{plot.regcorr}} for plotting the fitted correlations.
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
#' fitted(fit)
#'
#' @importFrom stats fitted
#' @export
fitted.regcorr <- function(object, ...) {
  return(object$fitted.rho)
}
