.fit_regcorr_engine <- function(Y, X, type, init, link, control,
                                warn = FALSE) {
  if (type == "normal") {
    NRfitBivNormal(
      Y = Y, X = X, betaIni = init, link = link, control = control
    )
  } else {
    NRfitBivBernoulli(
      Y = Y, X = X, beta0 = init, link = link, warn = warn,
      control = control
    )
  }
}

.bootstrap_fit_is_valid <- function(fit, type, p, control) {
  if (!is.list(fit) || !isTRUE(fit$converged) ||
      length(fit$betaCurrent) != p || any(!is.finite(fit$betaCurrent)) ||
      length(fit$logLik) != 1L || !is.finite(fit$logLik)) {
    return(FALSE)
  }

  if (type == "normal") {
    return(length(fit$rho) > 0L && all(is.finite(fit$rho)) &&
             all(abs(fit$rho) < 1 - control$boundary_eps))
  }

  length(fit$minJointProbability) == 1L &&
    is.finite(fit$minJointProbability) &&
    fit$minJointProbability > control$boundary_eps
}

.bootstrap_regcorr <- function(Y, X, nboot, type, init, link, control,
                               fit_fun = .fit_regcorr_engine) {
  p <- ncol(X)
  coefficient_names <- colnames(X)
  vcov_matrix <- matrix(
    NA_real_, p, p,
    dimnames = list(coefficient_names, coefficient_names)
  )
  if (nboot == 0L) {
    return(list(
      vcov = vcov_matrix,
      nboot.valid = 0L,
      nboot.failed = 0L,
      nboot.nonconverged = 0L,
      nboot.errors = 0L,
      diagnostics = list(status = character(), error.messages = character())
    ))
  }

  n <- nrow(Y)
  boot_betas <- matrix(NA_real_, nboot, p)
  status <- rep("nonconverged", nboot)
  error_messages <- rep(NA_character_, nboot)

  for (i in seq_len(nboot)) {
    boot_fit <- tryCatch({
      index <- sample.int(n, n, replace = TRUE)
      fit_fun(
        Y = Y[index, , drop = FALSE],
        X = X[index, , drop = FALSE],
        type = type,
        init = init,
        link = link,
        control = control,
        warn = FALSE
      )
    }, error = function(e) e)

    if (inherits(boot_fit, "error")) {
      status[i] <- "error"
      error_messages[i] <- conditionMessage(boot_fit)
    } else if (.bootstrap_fit_is_valid(boot_fit, type, p, control)) {
      status[i] <- "valid"
      boot_betas[i, ] <- as.vector(boot_fit$betaCurrent)
    }
  }

  nboot_valid <- sum(status == "valid")
  nboot_nonconverged <- sum(status == "nonconverged")
  nboot_errors <- sum(status == "error")
  nboot_failed <- nboot_nonconverged + nboot_errors

  if (nboot_valid > p + 1L) {
    vcov_matrix <- stats::cov(
      boot_betas[status == "valid", , drop = FALSE]
    )
    dimnames(vcov_matrix) <- list(coefficient_names, coefficient_names)
  }

  failure_detail <- paste0(
    nboot_nonconverged, " optimizer non-convergence",
    if (nboot_nonconverged != 1L) "s" else "",
    "; ", nboot_errors, " error", if (nboot_errors != 1L) "s" else ""
  )
  if (nboot_valid <= p + 1L) {
    warning(
      "Only ", nboot_valid, " of ", nboot,
      " bootstrap replications were valid (", failure_detail,
      "); too few to estimate the variance-covariance matrix (set to NA).",
      call. = FALSE
    )
  } else if (nboot_failed > 0L) {
    warning(
      nboot_failed, " of ", nboot,
      " bootstrap replications were discarded (", failure_detail, ").",
      call. = FALSE
    )
  }

  names(status) <- paste0("replicate", seq_len(nboot))
  names(error_messages) <- names(status)
  error_messages <- error_messages[!is.na(error_messages)]
  list(
    vcov = vcov_matrix,
    nboot.valid = as.integer(nboot_valid),
    nboot.failed = as.integer(nboot_failed),
    nboot.nonconverged = as.integer(nboot_nonconverged),
    nboot.errors = as.integer(nboot_errors),
    diagnostics = list(
      status = status,
      error.messages = error_messages
    )
  )
}
