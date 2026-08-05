.correlation_link <- function(eta, link) {
  if (link == "1") stats::plogis(eta) else tanh(eta)
}

.invalid_objective <- function(message) {
  list(valid = FALSE, message = message)
}

.normal_marginal_statistics <- function(Y, X) {
  marginal_fit <- tryCatch(
    stats::lm.fit(X, Y),
    error = function(e) e
  )
  if (inherits(marginal_fit, "error")) {
    return(list(
      ok = FALSE,
      message = paste0("Normal marginal fitting failed: ",
                       conditionMessage(marginal_fit))
    ))
  }

  residual <- marginal_fit$residuals
  sigma_hat <- sqrt(colMeans(residual^2))
  scale_tol <- sqrt(.Machine$double.eps)
  if (any(!is.finite(residual)) || any(!is.finite(sigma_hat)) ||
      any(sigma_hat <= scale_tol)) {
    return(list(
      ok = FALSE,
      message = paste0(
        "Normal marginal residual scales are zero or non-finite; ",
        "the correlation likelihood is not estimable."
      )
    ))
  }

  y_tilde <- sweep(residual, 2L, sigma_hat, "/")
  list(
    ok = TRUE,
    sigma = sigma_hat,
    T1 = rowSums(y_tilde^2),
    T2 = y_tilde[, 1L] * y_tilde[, 2L]
  )
}

#' Score and Hessian weights for the bivariate normal likelihood
#'
#' Returns the observation-level score and Hessian weights with respect to the
#' linear predictor after the marginal estimates have been fixed. The formulas
#' are those in Dufera, Liu, and Xu (2023), Table 1.
#'
#' @param T1,T2 Fixed observation-level normal sufficient statistics.
#' @param rho Link-transformed correlations.
#' @param link Link function indicator (logistic or tanh).
#' @noRd
.normal_derivatives <- function(T1, T2, rho, link) {
  if (link == "1") {
    a <- -rho^2 / ((1 + rho)^2 * (1 - rho))
    b <- rho * (1 + rho^2) / ((1 + rho)^2 * (1 - rho))
    c <- rho^2 / (1 + rho)
    u <- -rho^2 * (2 - rho + rho^2) / ((1 + rho)^3 * (1 - rho))
    v <- -rho * (-1 + rho - 5 * rho^2 + rho^3) /
      ((1 + rho)^3 * (1 - rho))
    w <- rho^2 * (2 + rho) * (1 - rho) / (1 + rho)^2
  } else {
    a <- -rho / ((1 - rho) * (1 + rho))
    b <- (1 + rho^2) / ((1 - rho) * (1 + rho))
    c <- rho
    u <- -(1 + rho^2) / ((1 - rho) * (1 + rho))
    v <- 4 * rho / ((1 - rho) * (1 + rho))
    w <- (1 - rho) * (1 + rho)
  }

  list(
    score_weight = a * T1 + b * T2 + c,
    hessian_weight = u * T1 + v * T2 + w
  )
}

.normal_loglik_contributions <- function(T1, T2, rho) {
  -0.5 * log1p(-rho^2) -
    (T1 - 2 * rho * T2) / (2 * (1 - rho^2))
}

.normal_loglik <- function(T1, T2, rho) {
  sum(.normal_loglik_contributions(T1, T2, rho))
}

.normal_objective <- function(beta, X, T1, T2, link,
                              boundary_eps = 1e-10) {
  if (length(beta) != ncol(X) || any(!is.finite(beta))) {
    return(.invalid_objective("The coefficient vector is non-finite."))
  }

  rho <- as.vector(.correlation_link(X %*% beta, link))
  if (any(!is.finite(rho)) || any(abs(rho) >= 1 - boundary_eps)) {
    return(.invalid_objective(
      "At least one fitted normal correlation is on the numerical boundary."
    ))
  }

  loglik <- .normal_loglik(T1, T2, rho)
  derivatives <- .normal_derivatives(T1, T2, rho, link)
  score <- as.vector(crossprod(X, derivatives$score_weight))
  hessian <- crossprod(
    X,
    X * as.vector(derivatives$hessian_weight)
  )

  if (!is.finite(loglik) || any(!is.finite(score)) ||
      any(!is.finite(hessian))) {
    return(.invalid_objective(
      "The normal log-likelihood or its derivatives are non-finite."
    ))
  }

  list(
    valid = TRUE,
    loglik = loglik,
    score = score,
    hessian = hessian,
    rho = rho,
    min_joint_probability = NA_real_
  )
}

#' Score and Hessian weights for one cell of the bivariate Bernoulli likelihood
#'
#' For a cell with log-probability log(c + d*rho), rho = h(eta), returns the
#' score weight e = d*h'(eta)/(c + d*rho) and the Hessian weight
#' f = -d^2*(h'(eta))^2/(c + d*rho)^2 + d*h''(eta)/(c + d*rho).
#'
#' @param c Cell probability offset (vector).
#' @param d Cell probability slope in rho (vector).
#' @param rho Link-transformed correlation h(eta) (vector).
#' @param link Link function indicator (logistic or tanh).
#' @noRd
.bernoulli_derivatives <- function(c, d, rho, link) {
  if (link == "1") {
    hp <- rho * (1 - rho)
    hpp <- rho * (1 - rho) * (1 - 2 * rho)
  } else {
    hp <- 1 - rho^2
    hpp <- -2 * rho * (1 - rho^2)
  }
  denom <- c + d * rho
  list(
    score_weight = d * hp / denom,
    hessian_weight = -d^2 * hp^2 / denom^2 + d * hpp / denom
  )
}

.bernoulli_likelihood_data <- function(p1, p2, Y) {
  p1 <- as.vector(p1)
  p2 <- as.vector(p2)
  d11 <- sqrt(p1 * (1 - p1) * p2 * (1 - p2))

  list(
    c00 = (1 - p1) * (1 - p2),
    c01 = (1 - p1) * p2,
    c10 = p1 * (1 - p2),
    c11 = p1 * p2,
    d11 = d11,
    I00 = (1 - Y[, 1L]) * (1 - Y[, 2L]),
    I01 = (1 - Y[, 1L]) * Y[, 2L],
    I10 = Y[, 1L] * (1 - Y[, 2L]),
    I11 = Y[, 1L] * Y[, 2L]
  )
}

.bernoulli_joint_probabilities <- function(likelihood_data, rho) {
  with(likelihood_data, cbind(
    p00 = c00 + d11 * rho,
    p01 = c01 - d11 * rho,
    p10 = c10 - d11 * rho,
    p11 = c11 + d11 * rho
  ))
}

.bernoulli_feasible <- function(likelihood_data, rho, eps = 1e-10) {
  if (any(!is.finite(rho)) || any(abs(rho) >= 1 - eps)) {
    return(FALSE)
  }
  joint <- .bernoulli_joint_probabilities(likelihood_data, rho)
  all(is.finite(joint)) && all(joint > eps)
}

.bernoulli_objective <- function(beta, X, likelihood_data, link,
                                 boundary_eps = 1e-10) {
  if (length(beta) != ncol(X) || any(!is.finite(beta))) {
    return(.invalid_objective("The coefficient vector is non-finite."))
  }

  rho <- as.vector(.correlation_link(X %*% beta, link))
  if (!.bernoulli_feasible(likelihood_data, rho, boundary_eps)) {
    return(.invalid_objective(
      "At least one bivariate Bernoulli joint probability is infeasible."
    ))
  }

  joint <- .bernoulli_joint_probabilities(likelihood_data, rho)
  with(likelihood_data, {
    loglik <- sum(
      I00 * log(joint[, "p00"]) + I01 * log(joint[, "p01"]) +
        I10 * log(joint[, "p10"]) + I11 * log(joint[, "p11"])
    )

    a00 <- .bernoulli_derivatives(c00, d11, rho, link)
    a01 <- .bernoulli_derivatives(c01, -d11, rho, link)
    a10 <- .bernoulli_derivatives(c10, -d11, rho, link)
    a11 <- .bernoulli_derivatives(c11, d11, rho, link)

    score_weight <- I00 * a00$score_weight + I01 * a01$score_weight +
      I10 * a10$score_weight + I11 * a11$score_weight
    hessian_weight <- I00 * a00$hessian_weight +
      I01 * a01$hessian_weight + I10 * a10$hessian_weight +
      I11 * a11$hessian_weight

    score <- as.vector(crossprod(X, score_weight))
    hessian <- crossprod(X, X * as.vector(hessian_weight))

    if (!is.finite(loglik) || any(!is.finite(score)) ||
        any(!is.finite(hessian))) {
      return(.invalid_objective(
        "The Bernoulli log-likelihood or its derivatives are non-finite."
      ))
    }

    list(
      valid = TRUE,
      loglik = loglik,
      score = score,
      hessian = hessian,
      rho = rho,
      joint_probabilities = joint,
      min_joint_probability = min(joint)
    )
  })
}

.fit_bernoulli_margins <- function(Y, X, warn = FALSE) {
  fits <- tryCatch(
    list(
      fit1 = suppressWarnings(
        stats::glm.fit(X, Y[, 1L], family = stats::binomial(link = "logit"))
      ),
      fit2 = suppressWarnings(
        stats::glm.fit(X, Y[, 2L], family = stats::binomial(link = "logit"))
      )
    ),
    error = function(e) e
  )

  if (inherits(fits, "error")) {
    message <- paste0("Marginal logistic fitting failed: ",
                      conditionMessage(fits))
    if (warn) warning(message, call. = FALSE)
    return(list(ok = FALSE, message = message))
  }

  p1 <- as.vector(fits$fit1$fitted.values)
  p2 <- as.vector(fits$fit2$fitted.values)
  separated <- !isTRUE(fits$fit1$converged) || !isTRUE(fits$fit2$converged) ||
    any(p1 <= 1e-8) || any(p1 >= 1 - 1e-8) ||
    any(p2 <= 1e-8) || any(p2 >= 1 - 1e-8)

  if (warn && separated) {
    warning(
      "Marginal logistic fits produced fitted probabilities near 0 or 1 ",
      "or did not converge (perfect or quasi-perfect separation); the ",
      "correlation parameter may be weakly identified.",
      call. = FALSE
    )
  }

  if (length(p1) != nrow(Y) || length(p2) != nrow(Y) ||
      any(!is.finite(p1)) || any(!is.finite(p2)) ||
      any(p1 <= 0) || any(p1 >= 1) || any(p2 <= 0) || any(p2 >= 1)) {
    return(list(
      ok = FALSE,
      message = paste0(
        "Marginal logistic fitted probabilities are non-finite or on the ",
        "boundary; the correlation likelihood is not estimable."
      ),
      separated = separated
    ))
  }

  list(ok = TRUE, p1 = p1, p2 = p2, separated = separated)
}

.safe_condition_number <- function(hessian) {
  if (!is.matrix(hessian) || any(!is.finite(hessian))) return(Inf)
  value <- tryCatch(
    kappa(hessian, exact = TRUE),
    error = function(e) Inf
  )
  if (length(value) != 1L || !is.finite(value)) Inf else as.numeric(value)
}

.optimizer_result <- function(state, beta, num_iter, converged, message,
                              step_size = NA_real_, relative_change = NA_real_,
                              condition_number = NA_real_, num_halving = 0L,
                              start_adjusted = FALSE) {
  score <- if (isTRUE(state$valid)) as.vector(state$score) else
    rep(NA_real_, length(beta))
  gradient_norm <- if (all(is.finite(score))) max(abs(score)) else Inf
  list(
    betaCurrent = as.vector(beta),
    numIter = as.integer(num_iter),
    restart = 0L,
    converged = isTRUE(converged),
    finalScore = score,
    gradientNorm = gradient_norm,
    stepSize = step_size,
    relativeChange = relative_change,
    convergenceMessage = message,
    logLik = if (isTRUE(state$valid)) state$loglik else NA_real_,
    conditionNumber = condition_number,
    numHalving = as.integer(num_halving),
    startAdjusted = isTRUE(start_adjusted),
    rho = if (isTRUE(state$valid)) state$rho else rep(NA_real_, 0L),
    minJointProbability = if (isTRUE(state$valid))
      state$min_joint_probability else NA_real_
  )
}

.failed_optimizer_result <- function(beta, message, start_adjusted = FALSE) {
  .optimizer_result(
    state = .invalid_objective(message),
    beta = beta,
    num_iter = 0L,
    converged = FALSE,
    message = message,
    start_adjusted = start_adjusted
  )
}

.find_admissible_start <- function(beta, anchor, evaluate, control) {
  initial_state <- tryCatch(
    evaluate(beta),
    error = function(e) .invalid_objective(conditionMessage(e))
  )
  if (isTRUE(initial_state$valid)) {
    return(list(ok = TRUE, beta = beta, state = initial_state,
                adjusted = FALSE))
  }

  anchor_state <- tryCatch(
    evaluate(anchor),
    error = function(e) .invalid_objective(conditionMessage(e))
  )
  if (!isTRUE(anchor_state$valid)) {
    return(list(
      ok = FALSE,
      beta = beta,
      state = initial_state,
      adjusted = FALSE,
      message = paste0(
        "Neither the supplied initial value nor a deterministic interior ",
        "anchor is admissible: ", anchor_state$message
      )
    ))
  }

  if (all(is.finite(beta))) {
    for (halving in seq_len(control$max_halving)) {
      scale <- 2^(-halving)
      if (scale < control$min_step) break
      candidate <- anchor + scale * (beta - anchor)
      candidate_state <- tryCatch(
        evaluate(candidate),
        error = function(e) .invalid_objective(conditionMessage(e))
      )
      if (isTRUE(candidate_state$valid)) {
        return(list(ok = TRUE, beta = candidate, state = candidate_state,
                    adjusted = TRUE))
      }
    }
  }

  list(ok = TRUE, beta = anchor, state = anchor_state, adjusted = TRUE)
}

.newton_optimize <- function(beta, evaluate, control, state = NULL,
                             start_adjusted = FALSE) {
  if (is.null(state)) {
    state <- tryCatch(
      evaluate(beta),
      error = function(e) .invalid_objective(conditionMessage(e))
    )
  }
  if (!isTRUE(state$valid)) {
    return(.failed_optimizer_result(
      beta,
      paste0("Initial objective evaluation failed: ", state$message),
      start_adjusted = start_adjusted
    ))
  }

  gradient_norm <- max(abs(state$score))
  if (gradient_norm < control$gradtol) {
    return(.optimizer_result(
      state, beta, 0L, TRUE,
      "The score tolerance was satisfied at the initial value.",
      step_size = 0, relative_change = 0,
      condition_number = .safe_condition_number(state$hessian),
      start_adjusted = start_adjusted
    ))
  }

  last_step_size <- NA_real_
  relative_change <- NA_real_
  condition_number <- NA_real_
  num_halving <- 0L

  for (iteration in seq_len(control$maxit)) {
    condition_number <- .safe_condition_number(state$hessian)
    if (condition_number > control$hessian_max_condition) {
      return(.optimizer_result(
        state, beta, iteration, FALSE,
        paste0(
          "The Hessian is singular or ill-conditioned (condition number ",
          format(condition_number, digits = 4), ")."
        ),
        step_size = last_step_size,
        relative_change = relative_change,
        condition_number = condition_number,
        num_halving = num_halving,
        start_adjusted = start_adjusted
      ))
    }

    newton_step <- tryCatch(
      solve(state$hessian, state$score),
      error = function(e) e
    )
    if (inherits(newton_step, "error") || any(!is.finite(newton_step))) {
      detail <- if (inherits(newton_step, "error"))
        conditionMessage(newton_step) else "the Newton direction is non-finite"
      return(.optimizer_result(
        state, beta, iteration, FALSE,
        paste0("The Newton system could not be solved: ", detail, "."),
        step_size = last_step_size,
        relative_change = relative_change,
        condition_number = condition_number,
        num_halving = num_halving,
        start_adjusted = start_adjusted
      ))
    }

    accepted <- FALSE
    candidate <- NULL
    candidate_state <- NULL
    halving_tried <- 0L
    objective_tol <- 100 * .Machine$double.eps *
      (1 + abs(state$loglik))

    for (halving in 0:control$max_halving) {
      step_size <- 2^(-halving)
      if (step_size < control$min_step) break
      halving_tried <- halving
      candidate <- as.vector(beta - step_size * newton_step)
      if (all(is.finite(candidate))) {
        candidate_state <- tryCatch(
          evaluate(candidate),
          error = function(e) .invalid_objective(conditionMessage(e))
        )
        accepted <- isTRUE(candidate_state$valid) &&
          candidate_state$loglik >= state$loglik - objective_tol
      }
      if (accepted) break
    }

    num_halving <- num_halving + halving_tried
    if (!accepted) {
      return(.optimizer_result(
        state, beta, iteration, FALSE,
        paste0(
          "Step-halving could not find an admissible, non-decreasing Newton ",
          "update within the configured limits."
        ),
        step_size = last_step_size,
        relative_change = relative_change,
        condition_number = condition_number,
        num_halving = num_halving,
        start_adjusted = start_adjusted
      ))
    }

    relative_change <- max(abs(candidate - beta)) /
      (1 + max(abs(beta)))
    beta <- candidate
    state <- candidate_state
    last_step_size <- step_size
    gradient_norm <- max(abs(state$score))

    if (relative_change < control$reltol &&
        gradient_norm < control$gradtol) {
      return(.optimizer_result(
        state, beta, iteration, TRUE,
        "Both the relative-change and score tolerances were satisfied.",
        step_size = last_step_size,
        relative_change = relative_change,
        condition_number = .safe_condition_number(state$hessian),
        num_halving = num_halving,
        start_adjusted = start_adjusted
      ))
    }

    if (relative_change <= .Machine$double.eps &&
        gradient_norm >= control$gradtol) {
      return(.optimizer_result(
        state, beta, iteration, FALSE,
        "The Newton iteration stagnated above the score tolerance.",
        step_size = last_step_size,
        relative_change = relative_change,
        condition_number = .safe_condition_number(state$hessian),
        num_halving = num_halving,
        start_adjusted = start_adjusted
      ))
    }
  }

  .optimizer_result(
    state, beta, control$maxit, FALSE,
    "The maximum number of Newton iterations was reached.",
    step_size = last_step_size,
    relative_change = relative_change,
    condition_number = .safe_condition_number(state$hessian),
    num_halving = num_halving,
    start_adjusted = start_adjusted
  )
}

.bernoulli_anchor <- function(X, likelihood_data, link, boundary_eps) {
  p <- ncol(X)
  d11 <- likelihood_data$d11
  if (any(!is.finite(d11)) || any(d11 <= 0)) {
    return(NULL)
  }

  lower <- pmax(
    (boundary_eps - likelihood_data$c00) / d11,
    (boundary_eps - likelihood_data$c11) / d11
  )
  upper <- pmin(
    (likelihood_data$c01 - boundary_eps) / likelihood_data$d11,
    (likelihood_data$c10 - boundary_eps) / likelihood_data$d11
  )
  link_lower <- if (link == "1") 0 else -1 + boundary_eps
  feasible_lower <- max(c(lower, link_lower))
  feasible_upper <- min(c(upper, 1 - boundary_eps))
  if (!is.finite(feasible_lower) || !is.finite(feasible_upper) ||
      feasible_lower >= feasible_upper) {
    return(NULL)
  }

  target_rho <- feasible_lower + (feasible_upper - feasible_lower) / 2
  target_eta <- if (link == "1") stats::qlogis(target_rho) else
    atanh(target_rho)
  if (!is.finite(target_eta)) return(NULL)

  constant_columns <- which(vapply(seq_len(p), function(j) {
    all(is.finite(X[, j])) && abs(X[1L, j]) > sqrt(.Machine$double.eps) &&
      max(abs(X[, j] - X[1L, j])) <= sqrt(.Machine$double.eps)
  }, logical(1)))

  anchor <- rep(0, p)
  if (length(constant_columns)) {
    j <- constant_columns[1L]
    anchor[j] <- target_eta / X[1L, j]
    return(anchor)
  }

  approximate <- tryCatch(
    stats::lm.fit(X, rep(target_eta, nrow(X)))$coefficients,
    error = function(e) NULL
  )
  if (is.null(approximate)) return(NULL)
  approximate[!is.finite(approximate)] <- 0
  as.vector(approximate)
}

#' Estimate beta for Bivariate Normal responses using Newton Raphson
#'
#' @param Y n by 2 matrix, paired responses.
#' @param X n by p matrix, covariate matrix including first column of ones.
#' @param betaIni Initial estimate of beta.
#' @param link Indicator of link function ("1" = logistic, "2" = tanh).
#' @param control Numerical controls from \code{\link{regcorr_control}}.
#' @return A list containing the coefficient estimate, convergence status,
#'   iteration diagnostics, final score norm, accepted step size, and message.
#' @importFrom stats lm.fit
NRfitBivNormal <- function(Y, X, betaIni, link,
                           control = regcorr_control()) {
  control <- .validate_regcorr_control(control)
  marginal <- .normal_marginal_statistics(Y, X)
  if (!isTRUE(marginal$ok)) {
    return(.failed_optimizer_result(betaIni, marginal$message))
  }

  evaluate <- function(beta) {
    .normal_objective(
      beta, X, marginal$T1, marginal$T2, link,
      boundary_eps = control$boundary_eps
    )
  }
  start <- .find_admissible_start(
    betaIni, rep(0, ncol(X)), evaluate, control
  )
  if (!isTRUE(start$ok)) {
    return(.failed_optimizer_result(betaIni, start$message))
  }

  .newton_optimize(
    start$beta, evaluate, control, state = start$state,
    start_adjusted = start$adjusted
  )
}

#' Estimate beta for Bivariate Bernoulli responses using Newton Raphson
#'
#' @param Y n by 2 matrix, paired responses.
#' @param X n by p matrix, covariate matrix including first column of ones.
#' @param beta0 Initial estimate of beta.
#' @param link Indicator of link function ("1" = logistic, "2" = tanh).
#' @param warn Whether to warn when marginal logistic fits exhibit separation.
#' @param control Numerical controls from \code{\link{regcorr_control}}.
#' @return A list containing the coefficient estimate, convergence status,
#'   iteration diagnostics, final score norm, accepted step size, and message.
#' @importFrom stats glm.fit binomial
NRfitBivBernoulli <- function(Y, X, beta0, link, warn = FALSE,
                              control = regcorr_control()) {
  control <- .validate_regcorr_control(control)
  marginal <- .fit_bernoulli_margins(Y, X, warn = warn)
  if (!isTRUE(marginal$ok)) {
    return(.failed_optimizer_result(beta0, marginal$message))
  }

  likelihood_data <- .bernoulli_likelihood_data(
    marginal$p1, marginal$p2, Y
  )
  evaluate <- function(beta) {
    .bernoulli_objective(
      beta, X, likelihood_data, link,
      boundary_eps = control$boundary_eps
    )
  }
  anchor <- .bernoulli_anchor(
    X, likelihood_data, link, control$boundary_eps
  )
  if (is.null(anchor)) {
    return(.failed_optimizer_result(
      beta0,
      paste0(
        "No numerically interior Bernoulli starting value could be found ",
        "for the fitted marginal probabilities."
      )
    ))
  }

  start <- .find_admissible_start(beta0, anchor, evaluate, control)
  if (!isTRUE(start$ok)) {
    return(.failed_optimizer_result(beta0, start$message))
  }

  .newton_optimize(
    start$beta, evaluate, control, state = start$state,
    start_adjusted = start$adjusted
  )
}
