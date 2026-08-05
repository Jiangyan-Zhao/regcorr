make_optimizer_normal_data <- function(n = 400, link = "1", seed = 1) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  X <- cbind(`(Intercept)` = 1, x = x)
  beta <- if (link == "1") c(-0.5, 0.35) else c(0.05, 0.25)
  rho <- .correlation_link(X %*% beta, link)
  z1 <- rnorm(n)
  Y <- cbind(
    y1 = 0.4 + 0.3 * x + z1,
    y2 = -0.2 - 0.4 * x + rho * z1 +
      sqrt(1 - rho^2) * rnorm(n)
  )
  list(X = X, Y = Y, beta = beta, rho = as.vector(rho))
}

test_that("regcorr_control validates strict deterministic controls", {
  control <- regcorr_control()
  expect_equal(control$maxit, 100L)
  expect_equal(control$reltol, 1e-8)
  expect_equal(control$gradtol, 1e-6)
  expect_equal(control$max_halving, 20L)
  expect_s3_class(control, "regcorr_control")

  expect_error(regcorr_control(maxit = 0), "maxit")
  expect_error(regcorr_control(reltol = 0), "reltol")
  expect_error(regcorr_control(min_step = 2), "min_step")
  expect_error(.validate_regcorr_control(list(unknown = 1)), "Unknown")
})

test_that("normal optimizer converges for both links from imperfect starts", {
  for (link in c("1", "2")) {
    dat <- make_optimizer_normal_data(
      n = 500, link = link, seed = if (link == "1") 101 else 102
    )
    fit <- NRfitBivNormal(
      dat$Y, dat$X, betaIni = c(2.5, -2), link = link
    )

    expect_true(isTRUE(fit$converged), info = fit$convergenceMessage)
    expect_true(all(is.finite(fit$betaCurrent)))
    expect_lt(fit$gradientNorm, regcorr_control()$gradtol)
    expect_lt(fit$relativeChange, regcorr_control()$reltol)
    expect_true(all(abs(fit$rho) < 1 - regcorr_control()$boundary_eps))
    expect_equal(fit$restart, 0L)
  }
})

test_that("large normal fits are deterministic without random fallback", {
  dat <- make_optimizer_normal_data(n = 1500, link = "2", seed = 211)
  fit1 <- NRfitBivNormal(dat$Y, dat$X, c(0.1, 0.1), "2")
  fit2 <- NRfitBivNormal(dat$Y, dat$X, c(0.1, 0.1), "2")

  expect_true(fit1$converged)
  expect_true(fit2$converged)
  expect_identical(fit1$betaCurrent, fit2$betaCurrent)
  expect_identical(fit1$numIter, fit2$numIter)
  expect_identical(fit1$numHalving, fit2$numHalving)
  expect_equal(fit1$restart, 0L)
})

test_that("Bernoulli optimizer converges for logistic and tanh links", {
  settings <- list(
    list(link = "1", beta = c(-1, 0.3), seed = 301),
    list(link = "2", beta = c(0.05, 0.15), seed = 302)
  )
  for (setting in settings) {
    set.seed(setting$seed)
    dat <- genDataBB(
      numSample = 900, p = 1, betaTrue = setting$beta,
      eta1True = c(0.1, 0.2), eta2True = c(-0.1, -0.2),
      link = setting$link
    )
    fit <- NRfitBivBernoulli(
      dat$Y, dat$X, beta0 = c(0.1, 0.1), link = setting$link
    )

    expect_true(isTRUE(fit$converged), info = fit$convergenceMessage)
    expect_lt(fit$gradientNorm, regcorr_control()$gradtol)
    expect_gt(fit$minJointProbability, regcorr_control()$boundary_eps)
    expect_equal(fit$restart, 0L)
  }
})

test_that("Bernoulli step-halving recovers an infeasible full Newton step", {
  x <- seq(-1, 1, length.out = 30)
  X <- cbind(`(Intercept)` = 1, x = x)
  p1 <- plogis(0.6702294782543967 - 0.23594030441165614 * x)
  p2 <- plogis(0.17980893614059423 - 0.4540140675194353 * x)
  Y <- matrix(c(
    1,0, 1,1, 0,1, 1,0, 1,0, 1,1, 0,1, 1,0, 0,0, 1,0,
    0,0, 0,0, 1,1, 0,1, 1,0, 0,0, 1,0, 1,0, 1,1, 1,1,
    1,1, 0,0, 0,1, 1,1, 1,1, 0,0, 0,1, 0,0, 0,0, 0,1
  ), ncol = 2, byrow = TRUE)
  likelihood_data <- .bernoulli_likelihood_data(p1, p2, Y)
  control <- regcorr_control()
  evaluate <- function(beta) {
    .bernoulli_objective(
      beta, X, likelihood_data, "1", control$boundary_eps
    )
  }
  beta0 <- c(-1.5318688124367439, -3.0301567777071714)
  initial <- evaluate(beta0)
  full_step <- solve(initial$hessian, initial$score)

  expect_true(initial$valid)
  expect_false(evaluate(beta0 - full_step)$valid)

  fit <- .newton_optimize(beta0, evaluate, control, state = initial)
  expect_true(isTRUE(fit$converged), info = fit$convergenceMessage)
  expect_gte(fit$numHalving, 1L)
  expect_gt(fit$logLik, initial$loglik)
  expect_gt(fit$minJointProbability, control$boundary_eps)
})

test_that("Bernoulli feasibility checker handles near-boundary probabilities", {
  Y <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  likelihood_data <- .bernoulli_likelihood_data(
    rep(0.2, 4), rep(0.8, 4), Y
  )
  expect_true(.bernoulli_feasible(likelihood_data, rep(0.24, 4)))
  expect_false(.bernoulli_feasible(likelihood_data, rep(0.25, 4)))
})

test_that("Bernoulli anchor finds an interior point away from rho zero", {
  X <- matrix(1, 4, 1)
  Y <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  settings <- list(
    list(link = "1", p1 = 1e-8, p2 = 1e-8),
    list(link = "2", p1 = 1e-8, p2 = 1 - 1e-8)
  )

  for (setting in settings) {
    likelihood_data <- .bernoulli_likelihood_data(
      rep(setting$p1, 4), rep(setting$p2, 4), Y
    )
    anchor <- .bernoulli_anchor(X, likelihood_data, setting$link, 1e-10)
    expect_true(all(is.finite(anchor)))
    expect_true(.bernoulli_objective(
      anchor, X, likelihood_data, setting$link, 1e-10
    )$valid)
  }
})

test_that("degenerate Bernoulli data fail gracefully and warn on separation", {
  n <- 60
  x <- seq(-1, 1, length.out = n)
  y1 <- rep(0, n)
  y2 <- rep(c(0, 1), length.out = n)
  warnings_seen <- character()

  fit <- withCallingHandlers(
    regcorr(cbind(y1, y2) ~ x, type = "binary", nboot = 0),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(fit$converged)
  expect_true(all(is.finite(coef(fit))))
  expect_false(is.finite(fit$loglik))
  expect_true(any(grepl("Marginal logistic", warnings_seen)))
  expect_true(any(grepl("did not converge", warnings_seen)))
})

test_that("a converged Bernoulli fit has valid fitted joint probabilities", {
  set.seed(808)
  dat <- genDataBB(
    700, 1, betaTrue = c(-1, 0.2),
    eta1True = c(0.1, 0.2), eta2True = c(-0.1, -0.1), link = "1"
  )
  data <- data.frame(y1 = dat$Y[, 1], y2 = dat$Y[, 2], x = dat$X[, 2])
  fit <- regcorr(cbind(y1, y2) ~ x, data = data, nboot = 0)

  expect_true(fit$converged)
  expect_true(is.finite(fit$loglik))
  expect_gt(fit$min.joint.probability, fit$control$boundary_eps)

  p1 <- glm.fit(fit$x, fit$y[, 1], family = binomial("logit"))$fitted.values
  p2 <- glm.fit(fit$x, fit$y[, 2], family = binomial("logit"))$fitted.values
  likelihood_data <- .bernoulli_likelihood_data(p1, p2, fit$y)
  expect_true(.bernoulli_feasible(
    likelihood_data, fitted(fit), fit$control$boundary_eps
  ))
})
