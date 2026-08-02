#' Fit a Regression Model for the Pearson Correlation Coefficient
#'
#' Fits a regression model for a covariate-dependent Pearson correlation
#' coefficient between two responses. Both bivariate normal and bivariate
#' binary responses are supported; the model type is inferred from the
#' response by default (a response taking only values 0 and 1 is treated
#' as bivariate binary, otherwise bivariate normal) and can be set
#' explicitly via \code{type}.
#'
#' @param formula A model formula of the form \code{cbind(y1, y2) ~ x1 + x2}.
#' @param data Optional data frame containing the variables in \code{formula}.
#' @param type Model type: \code{"auto"} (default) infers the type from the
#'   response; \code{"normal"} and \code{"binary"} force a type.
#' @param link Correlation link function: \code{"logistic"} (or \code{"1"})
#'   or \code{"tanh"} (or \code{"2"}). Default is \code{"logistic"}.
#' @param init Optional initial values for the regression coefficients.
#'   Defaults to a vector of 0.1s.
#' @param nboot Number of bootstrap replications used to estimate the
#'   variance-covariance matrix of the coefficients. Set to 0 to skip the
#'   bootstrap (the resulting \code{vcov} will contain \code{NA}).
#'
#' @return An object of class \code{"regcorr_normal"} or
#'   \code{"regcorr_binary"} (and \code{"regcorr"}), with S3 methods
#'   \code{print}, \code{summary}, \code{predict}, \code{coef}, \code{vcov},
#'   \code{fitted}, \code{nobs}, and \code{plot}.
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
#' fit
#' summary(fit)
#'
#' @importFrom stats model.frame model.response model.matrix
#' @export
regcorr <- function(formula, data = NULL,
                    type = c("auto", "normal", "binary"),
                    link = c("logistic", "tanh", "1", "2"),
                    init = NULL, nboot = 100) {
  call <- match.call()

  # --- link ---
  link_arg <- match.arg(link)
  link_code <- switch(link_arg,
                      "logistic" = "1", "1" = "1",
                      "tanh"     = "2", "2" = "2")

  # --- parse formula and data ---
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)

  if (is.null(Y) || NCOL(Y) != 2) {
    stop("The response must be a 2-column matrix, e.g., cbind(y1, y2) ~ x1 + x2",
         call. = FALSE)
  }

  # --- model type ---
  type_arg <- match.arg(type)
  if (type_arg == "auto") {
    type_arg <- if (all(Y %in% c(0, 1))) "binary" else "normal"
  }

  # --- initial values ---
  p <- ncol(X)
  if (is.null(init)) init <- rep(0.1, p)
  if (length(init) != p) {
    stop("`init` must have length ", p, ", but has length ", length(init),
         call. = FALSE)
  }

  # --- point estimate via Newton-Raphson ---
  if (type_arg == "normal") {
    fit <- NRfitBivNormal(Y = Y, X = X, betaIni = init, link = link_code)
  } else {
    fit <- NRfitBivBernoulli(Y = Y, X = X, beta0 = init, link = link_code)
  }
  coefficients <- as.vector(fit$betaCurrent)
  names(coefficients) <- colnames(X)

  # --- bootstrap variance ---
  vcov_mat <- matrix(NA, p, p)
  if (nboot > 0) {
    n <- nrow(Y)
    boot_betas <- matrix(0, nboot, p)
    for (i in seq_len(nboot)) {
      idx <- sample(n, replace = TRUE)
      if (type_arg == "normal") {
        boot_fit <- NRfitBivNormal(Y = Y[idx, , drop = FALSE],
                                   X = X[idx, , drop = FALSE],
                                   betaIni = init, link = link_code)
      } else {
        boot_fit <- NRfitBivBernoulli(Y = Y[idx, , drop = FALSE],
                                      X = X[idx, , drop = FALSE],
                                      beta0 = init, link = link_code)
      }
      boot_betas[i, ] <- as.vector(boot_fit$betaCurrent)
    }
    vcov_mat <- stats::cov(boot_betas)
    colnames(vcov_mat) <- rownames(vcov_mat) <- colnames(X)
  }

  # --- fitted correlations ---
  rho_fitted <- switch(link_code,
                       "1" = logistic(X %*% coefficients),
                       "2" = tanh(X %*% coefficients))

  res <- list(
    coefficients = coefficients,
    vcov         = vcov_mat,
    fitted.rho   = as.vector(rho_fitted),
    type         = type_arg,
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
  class(res) <- c(paste0("regcorr_", type_arg), "regcorr")
  res
}
