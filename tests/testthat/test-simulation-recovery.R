test_that("normal simulation recovers the correlation trend", {
  set.seed(1201)
  beta_true <- c(0.1, 0.35)
  dat <- genDataBN(
    numSample = 1000, p = 1, betaTrue = beta_true,
    eta1True = c(0, 0), eta2True = c(0, 0), link = "2"
  )
  frame <- data.frame(
    y1 = dat$Y[, 1], y2 = dat$Y[, 2], x = dat$X[, 2]
  )
  fit <- regcorr(
    cbind(y1, y2) ~ x, data = frame, type = "normal",
    link = "tanh", nboot = 0
  )

  expect_true(fit$converged)
  expect_true(all(is.finite(coef(fit))))
  expect_equal(unname(coef(fit)), beta_true, tolerance = 0.45)
  expect_gt(coef(fit)[["x"]], 0)
  expect_true(all(abs(fitted(fit)) < 1))
})

test_that("Bernoulli simulation recovers the positive correlation trend", {
  set.seed(1202)
  beta_true <- c(-1, 0.5)
  dat <- genDataBB(
    numSample = 1800, p = 1, betaTrue = beta_true,
    eta1True = c(0, 0), eta2True = c(0, 0), link = "1"
  )
  frame <- data.frame(
    y1 = dat$Y[, 1], y2 = dat$Y[, 2], x = dat$X[, 2]
  )
  fit <- regcorr(
    cbind(y1, y2) ~ x, data = frame, type = "binary",
    link = "logistic", nboot = 0
  )

  expect_true(fit$converged)
  expect_true(all(is.finite(coef(fit))))
  expect_equal(unname(coef(fit)), beta_true, tolerance = 0.65)
  expect_gt(coef(fit)[["x"]], 0)
  expect_true(all(fitted(fit) > 0 & fitted(fit) < 1))
  expect_gt(fit$min.joint.probability, fit$control$boundary_eps)
})
