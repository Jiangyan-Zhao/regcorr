#' Print a Fitted Correlation Regression Model
#'
#' Prints a fitted \code{regcorr} model: the original call, the estimated
#' correlation-link coefficients, and optimization details (number of
#' Newton-Raphson iterations and restarts).
#'
#' @param x An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The argument \code{x}, returned invisibly.
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{summary.regcorr}} for the coefficient table with bootstrap
#'   standard errors; \code{\link{coef.regcorr}} for extracting the
#'   coefficients.
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
#' print(fit)
#'
#' @export
print.regcorr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients (Correlation Link Beta):\n")
  print(x$coefficients)
  cat("\nOptimization: Iterations =", x$numIter, "| Restarts =", x$restart, "\n")
  invisible(x)
}

#' Print a Summary of a Correlation Regression Model
#'
#' Prints a \code{summary.regcorr} object produced by
#' \code{\link{summary.regcorr}}, including the model type, link function,
#' bootstrap standard errors, and the coefficient table with z-statistics
#' and p-values.
#'
#' @param x An object of class \code{"summary.regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The argument \code{x}, returned invisibly.
#'
#' @seealso \code{\link{summary.regcorr}} to create the summary;
#'   \code{\link{print.regcorr}} to print the fitted model itself.
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
#' print(summary(fit))
#'
#' @export
print.summary.regcorr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Model Type   :", x$class, "\n")
  cat("Link Function:", x$link, "\n")
  cat("Bootstrap S.E. (nboot = ", x$nboot, ")\n\n", sep = "")

  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n---")
  cat("\nConvergence info: Iterations =", x$numIter, "| Restarts =", x$restart, "\n")
  invisible(x)
}
