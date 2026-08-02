#' Predict Correlations for New Observations
#'
#' Computes fitted or predicted Pearson correlation coefficients from a
#' fitted \code{regcorr} model. With \code{newdata = NULL} the fitted
#' correlations for the training data are returned; otherwise the
#' correlations are computed for the covariate values in \code{newdata}
#' using the inverse of the model's link function.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param newdata Optional data frame containing the covariate values for
#'   which correlations should be predicted. If \code{NULL}, the fitted
#'   correlations from the original data are returned.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predicted correlation coefficients.
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{fitted.regcorr}} for the fitted correlations of the
#'   training data; \code{\link{coef.regcorr}} for the coefficients on the
#'   link scale.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' x <- runif(n)
#' rho <- plogis(0.5 + 0.2 * x)
#' z1 <- rnorm(n)
#' y1 <- z1
#' y2 <- z1 * rho + rnorm(n) * sqrt(1 - rho^2)
#' dat <- data.frame(y1 = y1, y2 = y2, x = x)
#' fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 20)
#'
#' # Fitted correlations (same as fitted(fit))
#' predict(fit)
#'
#' # Predict at new covariate values
#' predict(fit, newdata = data.frame(x = c(0.1, 0.5, 0.9)))
#'
#' @importFrom stats delete.response model.frame model.matrix predict
#' @export
predict.regcorr <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.rho)
  }

  tt <- delete.response(object$terms)
  mf <- model.frame(tt, newdata)
  X_new <- model.matrix(tt, mf)

  link_code <- if (object$link %in% c("logistic", "1")) "1" else "2"

  pred_rho <- switch(link_code,
                     "1" = logistic(X_new %*% object$coefficients),
                     "2" = tanh(X_new %*% object$coefficients))
  return(as.vector(pred_rho))
}
