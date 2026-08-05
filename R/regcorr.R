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
#' @param control Numerical controls created by
#'   \code{\link{regcorr_control}}. A named
#'   list containing a subset of those controls is also accepted.
#'
#' @return An object of class \code{"regcorr_normal"} or
#'   \code{"regcorr_binary"} (and \code{"regcorr"}), with S3 methods
#'   \code{print}, \code{summary}, \code{predict}, \code{coef}, \code{vcov},
#'   \code{fitted}, \code{nobs}, and \code{plot}. The object also contains
#'   \itemize{
#'     \item \code{converged}: whether both the relative coefficient-change
#'       and score tolerances were met;
#'     \item \code{point.objective.valid}: whether the optimizer established a
#'       finite, model-domain-valid final objective state;
#'     \item \code{margins}: fitted nuisance-model quantities retained for
#'       inference. Normal fits contain marginal coefficients, residual scales,
#'       and residuals; binary fits contain both marginal-logistic coefficient
#'       vectors, fitted probabilities, and convergence/separation diagnostics;
#'     \item \code{gradient.norm}, \code{step.size}, and
#'       \code{convergence.message}: final optimizer diagnostics;
#'     \item \code{xlevels} and \code{contrasts}: training factor levels and
#'       contrast specifications used to construct compatible prediction design
#'       matrices;
#'     \item \code{na.action}: the \code{na.action} used (or \code{NULL});
#'     \item \code{nboot.valid}: the number of bootstrap replications
#'       retained after discarding non-converged and errored fits;
#'     \item \code{nboot.failed}, \code{nboot.nonconverged}, and
#'       \code{nboot.errors}: bootstrap failure diagnostics.
#'     \item \code{bootstrap.skipped} and \code{bootstrap.skip.reason}:
#'       whether requested bootstrap inference was skipped because the point
#'       objective was invalid, and the corresponding reason.
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
                    na.action = getOption("na.action"),
                    control = regcorr_control()) {
  call <- match.call()
  control <- .validate_regcorr_control(control)

  if (!is.numeric(nboot) || length(nboot) != 1L || is.na(nboot) ||
      !is.finite(nboot) || nboot < 0 || nboot > .Machine$integer.max ||
      nboot != floor(nboot)) {
    stop("`nboot` must be one non-negative whole number.", call. = FALSE)
  }
  nboot <- as.integer(nboot)

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
  xlevels <- stats::.getXlevels(mt, mf)
  contrasts <- attr(X, "contrasts")

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
  if (!is.numeric(Y) || any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("The response and covariates must contain only finite numeric values.",
         call. = FALSE)
  }

  # --- model type ---
  type_arg <- match.arg(type)
  if (type_arg == "auto") {
    type_arg <- if (all(Y %in% c(0, 1))) "binary" else "normal"
  }
  if (type_arg == "binary" && !all(Y %in% c(0, 1))) {
    stop("For `type = \"binary\"`, both response columns must contain only 0 and 1.",
         call. = FALSE)
  }

  # --- initial values ---
  p <- ncol(X)
  if (is.null(init)) init <- rep(0.1, p)
  if (length(init) != p) {
    stop("`init` must have length ", p, ", but has length ", length(init),
         call. = FALSE)
  }
  if (!is.numeric(init) || any(!is.finite(init))) {
    stop("`init` must contain only finite numeric values.", call. = FALSE)
  }
  init <- as.vector(init)

  # --- point estimate via safeguarded Newton-Raphson ---
  fit <- .fit_regcorr_engine(
    Y = Y, X = X, type = type_arg, init = init, link = link_code,
    control = control, warn = TRUE
  )
  coefficients <- as.vector(fit$betaCurrent)
  names(coefficients) <- colnames(X)
  final_score <- as.vector(fit$finalScore)
  names(final_score) <- colnames(X)
  point_objective_valid <- .fit_has_valid_objective_state(
    fit, type = type_arg, n = nrow(X), p = p, control = control
  )

  if (!point_objective_valid) {
    warning(
      "Point estimator did not establish a valid final objective state",
      if (nboot > 0L) "; bootstrap inference was skipped" else "",
      ". ", fit$convergenceMessage,
      call. = FALSE
    )
  } else if (!isTRUE(fit$converged)) {
    warning("Safeguarded Newton estimation did not converge after ",
            fit$numIter, " iteration(s): ", fit$convergenceMessage,
            " The reported coefficients may be unreliable.",
            call. = FALSE)
  }

  # --- bootstrap variance ---
  bootstrap <- if (point_objective_valid) {
    .bootstrap_regcorr(
      Y = Y, X = X, nboot = nboot, type = type_arg, init = init,
      link = link_code, control = control
    )
  } else {
    .empty_bootstrap_result(
      X = X,
      nboot = nboot,
      skipped = nboot > 0L,
      reason = paste0(
        "The point estimator did not establish a valid final objective state: ",
        fit$convergenceMessage
      )
    )
  }

  # --- fitted correlations ---
  rho_fitted <- if (point_objective_valid) {
    as.vector(fit$rho)
  } else {
    rep(NA_real_, nrow(X))
  }

  res <- list(
    coefficients = coefficients,
    vcov         = bootstrap$vcov,
    fitted.rho   = as.vector(rho_fitted),
    type         = type_arg,
    point.objective.valid = point_objective_valid,
    margins      = fit$margins,
    numIter      = fit$numIter,
    restart      = fit$restart,
    converged    = fit$converged,
    score        = final_score,
    gradient.norm = fit$gradientNorm,
    step.size    = fit$stepSize,
    relative.change = fit$relativeChange,
    convergence.message = fit$convergenceMessage,
    loglik       = fit$logLik,
    condition.number = fit$conditionNumber,
    num.halving  = fit$numHalving,
    start.adjusted = fit$startAdjusted,
    min.joint.probability = fit$minJointProbability,
    link         = link_arg,
    nboot        = nboot,
    nboot.valid  = bootstrap$nboot.valid,
    nboot.failed = bootstrap$nboot.failed,
    nboot.nonconverged = bootstrap$nboot.nonconverged,
    nboot.errors = bootstrap$nboot.errors,
    bootstrap.diagnostics = bootstrap$diagnostics,
    bootstrap.skipped = isTRUE(bootstrap$diagnostics$skipped),
    bootstrap.skip.reason = bootstrap$diagnostics$skip.reason,
    control      = control,
    call         = call,
    terms        = mt,
    model        = mf,
    xlevels      = xlevels,
    contrasts    = contrasts,
    na.action    = na_action,
    x            = X,
    y            = Y
  )
  class(res) <- c(paste0("regcorr_", type_arg), "regcorr")
  res
}
