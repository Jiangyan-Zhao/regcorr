mock_normal_fit <- function(beta, n, converged = TRUE) {
  list(
    betaCurrent = beta,
    converged = converged,
    logLik = if (converged) -1 else NA_real_,
    rho = if (converged) rep(0.2, n) else numeric(),
    minJointProbability = NA_real_
  )
}

test_that("bootstrap errors are isolated and excluded from covariance", {
  X <- matrix(1, 12, 1, dimnames = list(NULL, "(Intercept)"))
  Y <- cbind(rnorm(12), rnorm(12))
  calls <- 0L
  fit_fun <- function(Y, X, type, init, link, control, warn) {
    calls <<- calls + 1L
    if (calls == 2L) stop("intentional replicate error")
    if (calls == 4L) return(mock_normal_fit(calls, nrow(Y), FALSE))
    mock_normal_fit(calls, nrow(Y), TRUE)
  }

  set.seed(41)
  expect_warning(
    bootstrap <- .bootstrap_regcorr(
      Y, X, nboot = 6L, type = "normal", init = 0,
      link = "1", control = regcorr_control(), fit_fun = fit_fun
    ),
    "discarded"
  )

  expect_equal(bootstrap$nboot.valid, 4L)
  expect_equal(bootstrap$nboot.failed, 2L)
  expect_equal(bootstrap$nboot.nonconverged, 1L)
  expect_equal(bootstrap$nboot.errors, 1L)
  expect_equal(as.numeric(bootstrap$vcov), var(c(1, 3, 5, 6)))
  expect_equal(sum(bootstrap$diagnostics$status == "error"), 1L)
  expect_equal(sum(bootstrap$diagnostics$status == "nonconverged"), 1L)
  expect_match(
    unname(bootstrap$diagnostics$error.messages),
    "intentional replicate error"
  )
})

test_that("too few valid bootstrap fits return NA covariance with warning", {
  X <- matrix(1, 10, 1, dimnames = list(NULL, "(Intercept)"))
  Y <- cbind(rnorm(10), rnorm(10))
  calls <- 0L
  fit_fun <- function(Y, X, type, init, link, control, warn) {
    calls <<- calls + 1L
    if (calls == 2L) stop("replicate failed")
    mock_normal_fit(calls, nrow(Y), converged = calls == 1L)
  }

  set.seed(42)
  expect_warning(
    bootstrap <- .bootstrap_regcorr(
      Y, X, nboot = 3L, type = "normal", init = 0,
      link = "1", control = regcorr_control(), fit_fun = fit_fun
    ),
    "too few"
  )

  expect_equal(bootstrap$nboot.valid, 1L)
  expect_equal(bootstrap$nboot.failed, 2L)
  expect_equal(bootstrap$nboot.nonconverged, 1L)
  expect_equal(bootstrap$nboot.errors, 1L)
  expect_true(all(is.na(bootstrap$vcov)))
})
