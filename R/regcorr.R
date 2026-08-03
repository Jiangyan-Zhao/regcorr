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
#' @param na.action A function (or quoted string) specifying how missing
#'   values in the data are handled; passed to
#'   \code{\link[stats]{model.frame}}. The default drops incomplete rows
#'   before fitting, and the dropped rows are recorded in the
#'   \code{na.action} component of the returned object. The fitting routine
#'   itself does not support missing values.
#'
#' @return An object of class \code{"regcorr_normal"} or
#'   \code{"regcorr_binary"} (and \code{"regcorr"}), with S3 methods
#'   \code{print}, \code{summary}, \code{predict}, \code{coef}, \code{vcov},
#'   \code{fitted}, \code{nobs}, and \code{plot}. The object also contains
#'   \itemize{
#'     \item \code{converged}: whether the Newton-Raphson iterations
#'       converged within the iteration/restart limits;
#'     \item \code{na.action}: the \code{na.action} used (or \code{NULL});
#'     \item \code{nboot.valid}: the number of bootstrap replications
#'       retained after discarding non-converged fits.
#'   }
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
                    init = NULL, nboot = 100,
                    na.action = getOption("na.action")) {
  call <- match.call()

  # --- link ---
  link_arg <- match.arg(link)
  link_code <- switch(link_arg,
                      "logistic" = "1", "1" = "1",
                      "tanh"     = "2", "2" = "2")

  # --- parse formula and data ---
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)
  na_action <- attr(mf, "na.action")

  if (is.null(Y) || NCOL(Y) != 2) {
    stop("The response must be a 2-column matrix, e.g., cbind(y1, y2) ~ x1 + x2",
         call. = FALSE)
  }
  if (anyNA(Y) || anyNA(X)) {
    stop("Missing values in the response or covariates are not supported by ",
         "the fitting routine; use the default na.action = na.omit to drop ",
         "incomplete rows.",
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
    fit <- NRfitBivBernoulli(Y = Y, X = X, beta0 = init, link = link_code,
                             warn = TRUE)
  }
  coefficients <- as.vector(fit$betaCurrent)
  names(coefficients) <- colnames(X)

  if (!isTRUE(fit$converged)) {
    warning("Newton-Raphson estimation did not converge after ", fit$numIter,
            " iteration(s) and ", fit$restart, " restart(s); the reported ",
            "coefficients may be unreliable.",
            call. = FALSE)
  }

  # --- bootstrap variance ---
  vcov_mat <- matrix(NA, p, p)
  nboot_valid <- 0L
  if (nboot > 0) {
    n <- nrow(Y)
    boot_betas <- matrix(0, nboot, p)
    boot_ok <- logical(nboot)
    for (i in seq_len(nboot)) {
      idx <- sample(n, replace = TRUE)
      if (type_arg == "normal") {
        boot_fit <- NRfitBivNormal(Y = Y[idx, , drop = FALSE],
                                   X = X[idx, , drop = FALSE],
                                   betaIni = init, link = link_code)
      } else {
        boot_fit <- NRfitBivBernoulli(Y = Y[idx, , drop = FALSE],
                                      X = X[idx, , drop = FALSE],
                                      beta0 = init, link = link_code,
                                      warn = FALSE)
      }
      boot_betas[i, ] <- as.vector(boot_fit$betaCurrent)
      # discard replications that did not converge or did not move from init
      boot_ok[i] <- isTRUE(boot_fit$converged) &&
        sum(abs(boot_betas[i, ] - init)) >= 0.01
    }
    nboot_valid <- sum(boot_ok)
    if (nboot_valid > p + 1) {
      vcov_mat <- stats::cov(boot_betas[boot_ok, , drop = FALSE])
      colnames(vcov_mat) <- rownames(vcov_mat) <- colnames(X)
    }
    if (nboot_valid < nboot) {
      if (nboot_valid <= p + 1) {
        warning("Only ", nboot_valid, " of ", nboot, " bootstrap ",
                "replications were usable; too few to estimate the ",
                "variance-covariance matrix (set to NA).",
                call. = FALSE)
      } else {
        warning(nboot - nboot_valid, " of ", nboot, " bootstrap ",
                "replications were discarded (non-convergence or no movement ",
                "from the starting values).",
                call. = FALSE)
      }
    }
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
    converged    = fit$converged,
    link         = link_arg,
    nboot        = nboot,
    nboot.valid  = nboot_valid,
    call         = call,
    terms        = mt,
    model        = mf,
    na.action    = na_action,
    x            = X,
    y            = Y
  )
  class(res) <- c(paste0("regcorr_", type_arg), "regcorr")
  res
}
