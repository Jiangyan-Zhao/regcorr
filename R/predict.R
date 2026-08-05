#' Predict Correlations for New Observations
#'
#' Computes fitted or predicted Pearson correlation coefficients from a
#' fitted \code{regcorr} model. With \code{newdata = NULL} the fitted
#' correlations for the training data are returned; otherwise the
#' correlations are computed for the covariate values in \code{newdata}
#' using the inverse of the model's link function. Prediction model matrices
#' reuse the factor levels and contrasts recorded during fitting.
#'
#' @param object An object of class \code{"regcorr"}.
#' @param newdata Optional data frame containing the covariate values for
#'   which correlations should be predicted. If \code{NULL}, the fitted
#'   correlations from the original data are returned.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector of predicted correlation coefficients, aligned
#'   with the rows of \code{newdata}; rows with missing covariate values
#'   return \code{NA}.
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
#' @name predict
#' @importFrom stats delete.response model.frame model.matrix predict
#' @export
predict.regcorr <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.rho)
  }

  tt <- delete.response(object$terms)
  # keep incomplete rows so predictions align with newdata; they become NA
  mf <- tryCatch(
    model.frame(
      tt, newdata, na.action = stats::na.pass, xlev = object$xlevels
    ),
    error = function(e) {
      stop(
        "Cannot construct the prediction model frame from `newdata`: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  X_new <- tryCatch(
    model.matrix(tt, mf, contrasts.arg = object$contrasts),
    error = function(e) {
      stop(
        "Cannot construct a training-compatible prediction design matrix: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  expected_columns <- names(object$coefficients)
  if (is.null(expected_columns)) expected_columns <- colnames(object$x)
  missing_columns <- setdiff(expected_columns, colnames(X_new))
  extra_columns <- setdiff(colnames(X_new), expected_columns)
  if (length(expected_columns) != ncol(X_new) || length(missing_columns) ||
      length(extra_columns)) {
    stop(
      "The prediction design matrix is incompatible with the training ",
      "design matrix",
      if (length(missing_columns)) paste0(
        "; missing column(s): ", paste(missing_columns, collapse = ", ")
      ) else "",
      if (length(extra_columns)) paste0(
        "; unexpected column(s): ", paste(extra_columns, collapse = ", ")
      ) else "",
      ".",
      call. = FALSE
    )
  }
  X_new <- X_new[, expected_columns, drop = FALSE]

  if (isFALSE(object$point.objective.valid)) {
    return(rep(NA_real_, nrow(X_new)))
  }

  link_code <- if (object$link %in% c("logistic", "1")) "1" else "2"

  pred_rho <- switch(link_code,
                     "1" = logistic(X_new %*% object$coefficients),
                     "2" = tanh(X_new %*% object$coefficients))
  pred_rho[!stats::complete.cases(X_new)] <- NA
  return(as.vector(pred_rho))
}
