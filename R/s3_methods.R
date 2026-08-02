#' @export
print.regcorr <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients (Correlation Link Beta):\n")
  print(x$coefficients)
  cat("\nOptimization: Iterations =", x$numIter, "| Restarts =", x$restart, "\n")
  invisible(x)
}

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
    link         = object$link,
    nboot        = object$nboot,
    class        = class(object)[1]
  )
  class(res) <- "summary.regcorr"
  return(res)
}

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

#' Predict correlation (rho) for new observations
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

#' @export
coef.regcorr <- function(object, ...) {
  return(object$coefficients)
}

#' @export
vcov.regcorr <- function(object, ...) {
  return(object$vcov)
}

#' @export
fitted.regcorr <- function(object, ...) {
  return(object$fitted.rho)
}

#' @export
nobs.regcorr <- function(object, ...) {
  return(nrow(object$y))
}

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