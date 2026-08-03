#' Summarize a Fitted Correlation Regression Model
#'
#' Computes a summary of a fitted \code{regcorr} model: a coefficient table
#' with bootstrap standard errors, z-statistics and two-sided p-values,
#' together with the model type, link function and convergence information.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.regcorr"} (a list) with
#'   components:
#'   \describe{
#'     \item{call}{the model call;}
#'     \item{coefficients}{a matrix with columns \code{Estimate},
#'       \code{Std. Error}, \code{z value} and \code{Pr(>|z|)};}
#'     \item{numIter}{number of Newton-Raphson iterations;}
#'     \item{restart}{number of restarts used during optimization;}
#'     \item{link}{the link function used;}
#'     \item{nboot}{number of bootstrap replications;}
#'     \item{nboot.valid}{number of bootstrap replications retained after
#'       discarding non-converged fits;}
#'     \item{converged}{whether the Newton-Raphson iterations converged;}
#'     \item{class}{the model class (\code{"regcorr_normal"} or
#'       \code{"regcorr_binary"}).}
#'   }
#'
#' @seealso \code{\link{regcorr}} for model fitting;
#'   \code{\link{print.summary.regcorr}} for printing the summary;
#'   \code{\link{vcov.regcorr}} and \code{\link{coef.regcorr}} for direct
#'   access to the variance-covariance matrix and the coefficients.
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
#' summary(fit)
#'
#' @name summary
#' @export
summary.regcorr <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z_val <- object$coefficients / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_val)))

  coef_matrix <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )

  res <- list(
    call         = object$call,
    coefficients = coef_matrix,
    numIter      = object$numIter,
    restart      = object$restart,
    converged    = object$converged,
    link         = object$link,
    nboot        = object$nboot,
    nboot.valid  = object$nboot.valid,
    class        = class(object)[1]
  )
  class(res) <- "summary.regcorr"
  return(res)
}
