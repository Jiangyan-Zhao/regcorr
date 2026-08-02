#' Plot Fitted Correlations
#'
#' Produces plots of the fitted correlations of a \code{regcorr} model.
#' With \code{which = 1} (default), the fitted correlations are plotted
#' against the observation index; with \code{which = 2}, the fitted
#' correlation curve is plotted against the first non-intercept covariate.
#'
#' @param x An object of class \code{"regcorr"}.
#' @param which Which plot to draw: \code{1} for fitted correlations against
#'   the observation index, or \code{2} for the correlation curve against the
#'   first covariate.
#' @param ... Additional arguments passed to the underlying \code{plot}
#'   function.
#'
#' @return The argument \code{x}, returned invisibly.
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{fitted.regcorr}} for the fitted correlations;
#'   \code{\link{predict.regcorr}} for predictions at new covariate values.
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
#' plot(fit)
#' plot(fit, which = 2)
#'
#' @importFrom graphics plot abline points
#' @importFrom grDevices adjustcolor
#' @export
plot.regcorr <- function(x, which = 1, ...) {
  rho <- x$fitted.rho
  if (which == 1) {
    plot(rho, type = "b", pch = 19, col = "steelblue",
         xlab = "Observation Index", ylab = "Fitted Correlation (rho)",
         main = paste("Fitted Correlations -", class(x)[1]), ...)
    abline(h = 0, lty = 2, col = "gray")
  } else if (which == 2) {
    if (ncol(x$x) > 1) {
      x_var <- x$x[, 2]
      ord <- order(x_var)
      plot(x_var[ord], rho[ord], type = "l", lwd = 2, col = "firebrick",
           xlab = colnames(x$x)[2], ylab = "Fitted Correlation (rho)",
           main = "Correlation Curve vs Covariate", ...)
      points(x_var, rho, pch = 16, col = adjustcolor("black", alpha.f = 0.4))
    } else {
      plot(rho, type = "p", xlab = "Index", ylab = "Fitted Correlation (rho)", ...)
    }
  }
  invisible(x)
}
