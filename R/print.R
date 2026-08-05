#' Print a Fitted Correlation Regression Model
#'
#' Prints a fitted \code{regcorr} model: the original call, the estimated
#' correlation-link coefficients, and optimization details (number of
#' safeguarded Newton iterations, step-halving, and final score norm). Also
#' prints a
#' \code{summary.regcorr} object: the model type, link function, bootstrap
#' standard errors, and the coefficient table with z-statistics and p-values.
#'
#' @param x An object of class \code{"regcorr"} or \code{"summary.regcorr"}.
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
#' print(summary(fit))
#'
#' @name print
#' @rdname print
#' @export
print.regcorr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients (Correlation Link Beta):\n")
  print(x$coefficients)
  cat("\nOptimization: Iterations =", x$numIter, "| Restarts =", x$restart,
      if (isFALSE(x$converged)) "| NOT CONVERGED", "\n")
  if (isFALSE(x$point.objective.valid)) {
    cat("Final objective state: INVALID\n")
  }
  cat("Final score norm =", format(x$gradient.norm, digits = 4),
      "| Last step size =", format(x$step.size, digits = 4),
      "| Step halvings =", x$num.halving, "\n")
  if (isFALSE(x$converged)) {
    cat("Reason:", x$convergence.message, "\n")
  }
  invisible(x)
}

#' @rdname print
#' @export
print.summary.regcorr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Model Type   :", x$class, "\n")
  cat("Link Function:", x$link, "\n")
  if (isTRUE(x$bootstrap.skipped)) {
    cat("Bootstrap S.E.: skipped\n")
    if (!is.null(x$bootstrap.skip.reason)) {
      cat("Reason:", x$bootstrap.skip.reason, "\n")
    }
    cat("\n")
  } else if (x$nboot > 0) {
    cat("Bootstrap S.E. (nboot = ", x$nboot, sep = "")
    if (!is.null(x$nboot.valid) && x$nboot.valid < x$nboot) {
      cat(", valid = ", x$nboot.valid,
          ", non-converged = ", x$nboot.nonconverged,
          ", errors = ", x$nboot.errors, ")", sep = "")
    } else {
      cat(")", sep = "")
    }
    cat("\n\n")
  } else {
    cat("Bootstrap S.E.: not computed (nboot = 0)\n\n")
  }

  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n---")
  cat("\nConvergence info: Iterations =", x$numIter, "| Restarts =", x$restart,
      "| Converged =", if (isTRUE(x$converged)) "TRUE" else "FALSE", "\n")
  cat("Final score norm =", format(x$gradient.norm, digits = 4),
      "| Last step size =", format(x$step.size, digits = 4),
      "| Step halvings =", x$num.halving, "\n")
  cat("Reason:", x$convergence.message, "\n")
  invisible(x)
}
