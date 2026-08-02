#' Bivariate Normal Correlation Regression Model
#'
#' Fits a correlation regression model for bivariate normal responses using S3 methods.
#'
#' @param formula Model formula, e.g., \code{cbind(y1, y2) ~ x1 + x2}.
#' @param data Optional data frame containing variables in formula.
#' @param link Link function for correlation: "logistic" (or "1") or "tanh" (or "2"). Default is "logistic".
#' @param betaIni Initial vector for beta coefficients. Default is NULL (zero vector).
#' @param nboot Number of bootstrap iterations for standard error estimation. Default is 100.
#' @return An object of class \code{"regcorr_normal"} and \code{"regcorr"}.
#' @export
normal <- function(formula, data = NULL, link = c("logistic", "tanh", "1", "2"), 
                   betaIni = NULL, nboot = 100) {
  call <- match.call()
  
  # 1. 解析 link 参数
  link_arg <- match.arg(link)
  link_code <- switch(link_arg,
                      "logistic" = "1", "1" = "1",
                      "tanh"     = "2", "2" = "2")
  
  # 2. 解析 formula 和 data
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)
  
  if (is.null(Y) || ncol(Y) != 2) {
    stop("The response variable must be a 2-column matrix, e.g., cbind(y1, y2) ~ x1 + x2")
  }
  
  p <- ncol(X)
  if (is.null(betaIni)) {
    betaIni <- rep(0.1, p)
  }
  
  # 3. 拟合点估计
  fit <- NRfitBivNormal(Y = Y, X = X, betaIni = betaIni, link = link_code)
  coefficients <- as.vector(fit$betaCurrent)
  names(coefficients) <- colnames(X)
  
  # 4. Bootstrap 计算方差协方差矩阵 (来自 subRoutineTest 逻辑)
  vcov_mat <- matrix(NA, p, p)
  if (nboot > 0) {
    nSample <- nrow(Y)
    boot_betas <- matrix(0, nboot, p)
    for (i in 1:nboot) {
      idx <- sample(nSample, replace = TRUE)
      boot_fit <- NRfitBivNormal(Y = Y[idx, ], X = X[idx, , drop = FALSE], 
                                 betaIni = betaIni, link = link_code)
      boot_betas[i, ] <- as.vector(boot_fit$betaCurrent)
    }
    vcov_mat <- stats::cov(boot_betas)
    colnames(vcov_mat) <- rownames(vcov_mat) <- colnames(X)
  }
  
  # 5. 计算拟合的 correlation (rho)
  rho_fitted <- switch(link_code,
                       "1" = logistic(X %*% coefficients),
                       "2" = tanh(X %*% coefficients))
  
  # 构建 S3 对象
  res <- list(
    coefficients = coefficients,
    vcov         = vcov_mat,
    fitted.rho   = as.vector(rho_fitted),
    numIter      = fit$numIter,
    restart      = fit$restart,
    link         = link_arg,
    nboot        = nboot,
    call         = call,
    terms        = mt,
    model        = mf,
    x            = X,
    y            = Y
  )
  
  class(res) <- c("regcorr_normal", "regcorr")
  return(res)
}