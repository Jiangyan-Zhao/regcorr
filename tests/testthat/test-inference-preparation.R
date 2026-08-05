make_prep_normal_data <- function(n = 400, seed = 1701) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  rho <- plogis(-0.4 + 0.35 * x)
  z1 <- rnorm(n)
  data.frame(
    y1 = 0.2 + 0.4 * x + z1,
    y2 = -0.1 - 0.2 * x + rho * z1 +
      sqrt(1 - rho^2) * rnorm(n),
    x = x
  )
}

test_that("invalid point objectives return NA fitted values and skip bootstrap", {
  n <- 60
  dat <- data.frame(
    y1 = rep(0, n),
    y2 = rep(c(0, 1), length.out = n),
    x = seq(-1, 1, length.out = n)
  )
  warnings_seen <- character()

  fit <- withCallingHandlers(
    regcorr(
      cbind(y1, y2) ~ x, data = dat, type = "binary", nboot = 5
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(fit$point.objective.valid)
  expect_false(is.finite(fit$loglik))
  expect_true(all(is.na(fitted(fit))))
  expect_true(all(is.na(vcov(fit))))
  expect_equal(fit$nboot.valid, 0L)
  expect_equal(fit$nboot.failed, 0L)
  expect_true(fit$bootstrap.skipped)
  expect_true(all(fit$bootstrap.diagnostics$status == "skipped"))
  expect_length(fit$bootstrap.diagnostics$status, 5L)
  expect_match(fit$bootstrap.skip.reason, "valid final objective state")
  expect_true(any(grepl("bootstrap inference was skipped", warnings_seen)))

  new_prediction <- predict(fit, newdata = data.frame(x = c(-0.5, 0.5)))
  expect_true(all(is.na(new_prediction)))
})

test_that("valid non-converged point states retain validated rho and bootstrap", {
  dat <- make_prep_normal_data(n = 300, seed = 1702)
  control <- regcorr_control(maxit = 1L, reltol = 1e-14, gradtol = 1e-14)
  warnings_seen <- character()

  fit <- withCallingHandlers(
    regcorr(
      cbind(y1, y2) ~ x, data = dat, nboot = 3L, control = control
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  engine_fit <- .fit_regcorr_engine(
    Y = fit$y, X = fit$x, type = "normal",
    init = rep(0.1, ncol(fit$x)), link = "1", control = control
  )

  expect_false(fit$converged)
  expect_true(fit$point.objective.valid)
  expect_true(is.finite(fit$loglik))
  expect_equal(fitted(fit), engine_fit$rho)
  expect_false(anyNA(fitted(fit)))
  expect_false(fit$bootstrap.skipped)
  expect_length(fit$bootstrap.diagnostics$status, 3L)
  expect_false(any(fit$bootstrap.diagnostics$status == "skipped"))
  expect_true(any(grepl("did not converge", warnings_seen)))
})

test_that("converged fits retain inverse-link fitted behavior", {
  dat <- make_prep_normal_data(n = 500, seed = 1703)
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 0)

  expect_true(fit$converged)
  expect_true(fit$point.objective.valid)
  expect_equal(
    fitted(fit),
    as.vector(plogis(fit$x %*% coef(fit)))
  )
})

test_that("normal marginal nuisance estimates are retained", {
  dat <- make_prep_normal_data(n = 350, seed = 1704)
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 0)

  expect_type(fit$margins, "list")
  expect_equal(dim(fit$margins$coefficients), c(ncol(fit$x), 2L))
  expect_equal(rownames(fit$margins$coefficients), colnames(fit$x))
  expect_equal(dim(fit$margins$residuals), dim(fit$y))
  expect_length(fit$margins$sigma, 2L)
  expect_true(all(is.finite(fit$margins$sigma)))
  expect_true(all(fit$margins$sigma > 0))
})

test_that("binary marginal nuisance estimates are retained", {
  set.seed(1705)
  generated <- genDataBB(
    numSample = 900, p = 1, betaTrue = c(-1, 0.25),
    eta1True = c(0.1, 0.2), eta2True = c(-0.1, -0.15), link = "1"
  )
  dat <- data.frame(
    y1 = generated$Y[, 1], y2 = generated$Y[, 2],
    x = generated$X[, 2]
  )
  fit <- regcorr(
    cbind(y1, y2) ~ x, data = dat, type = "binary", nboot = 0
  )

  expect_type(fit$margins, "list")
  expect_length(fit$margins$coefficients1, ncol(fit$x))
  expect_length(fit$margins$coefficients2, ncol(fit$x))
  expect_equal(names(fit$margins$coefficients1), colnames(fit$x))
  expect_equal(names(fit$margins$coefficients2), colnames(fit$x))
  expect_length(fit$margins$fitted1, nrow(fit$y))
  expect_length(fit$margins$fitted2, nrow(fit$y))
  expect_true(all(fit$margins$fitted1 > 0 & fit$margins$fitted1 < 1))
  expect_true(all(fit$margins$fitted2 > 0 & fit$margins$fitted2 < 1))
  expect_length(fit$margins$converged, 2L)
  expect_type(fit$margins$separated, "logical")
})

test_that("prediction reuses training factor levels and interaction contrasts", {
  set.seed(1706)
  n_per_level <- 220
  g <- factor(
    rep(c("A", "B", "C"), each = n_per_level),
    levels = c("A", "B", "C")
  )
  x <- runif(length(g), -1, 1)
  eta <- -0.5 + 0.2 * x + 0.15 * (g == "B") - 0.1 * (g == "C") +
    0.12 * x * (g == "B") - 0.08 * x * (g == "C")
  rho <- plogis(eta)
  z1 <- rnorm(length(g))
  dat <- data.frame(
    y1 = z1,
    y2 = rho * z1 + sqrt(1 - rho^2) * rnorm(length(g)),
    g = g,
    x = x
  )
  fit <- regcorr(cbind(y1, y2) ~ g * x, data = dat, nboot = 0)

  expect_equal(fit$xlevels$g, c("A", "B", "C"))
  expect_false(is.null(fit$contrasts$g))

  only_a <- data.frame(g = factor(c("A", "A")), x = c(-0.4, 0.6))
  prediction <- predict(fit, newdata = only_a)
  expect_length(prediction, 2L)
  expect_true(all(is.finite(prediction)))

  prediction_frame <- model.frame(
    delete.response(fit$terms), only_a,
    na.action = stats::na.pass, xlev = fit$xlevels
  )
  prediction_matrix <- model.matrix(
    delete.response(fit$terms), prediction_frame,
    contrasts.arg = fit$contrasts
  )
  expect_identical(colnames(prediction_matrix), names(coef(fit)))

  aligned <- predict(
    fit,
    newdata = data.frame(
      g = factor(c("A", "A", "A")),
      x = c(-0.5, NA, 0.5)
    )
  )
  expect_length(aligned, 3L)
  expect_true(is.na(aligned[2L]))
  expect_false(anyNA(aligned[c(1L, 3L)]))

  expect_error(
    predict(
      fit,
      newdata = data.frame(g = factor("D"), x = 0)
    ),
    "new level"
  )
})

test_that("a valid Bernoulli beta0 does not require a fallback anchor", {
  x <- seq(-1, 1, length.out = 40)
  X <- cbind(`(Intercept)` = 1, x = x)
  Y <- matrix(
    rep(c(0, 0, 0, 1, 1, 0, 1, 1), 10),
    ncol = 2, byrow = TRUE
  )
  likelihood_data <- .bernoulli_likelihood_data(
    rep(0.5, nrow(X)), rep(0.5, nrow(X)), Y
  )
  control <- regcorr_control(maxit = 1L)
  evaluate <- function(beta) {
    .bernoulli_objective(
      beta, X, likelihood_data, link = "1",
      boundary_eps = control$boundary_eps
    )
  }
  beta0 <- c(stats::qlogis(0.1), 0.05)
  anchor_calls <- 0L

  start <- .find_bernoulli_start(
    beta0,
    evaluate,
    anchor_fun = function() {
      anchor_calls <<- anchor_calls + 1L
      NULL
    },
    control = control
  )

  expect_true(start$ok)
  expect_false(start$adjusted)
  expect_equal(start$beta, beta0)
  expect_true(start$state$valid)
  expect_equal(anchor_calls, 0L)

  fit <- .newton_optimize(
    start$beta, evaluate, control, state = start$state,
    start_adjusted = start$adjusted
  )
  expect_true(is.finite(fit$logLik))
  expect_false(fit$startAdjusted)
})
