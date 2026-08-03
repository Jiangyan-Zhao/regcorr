# End-to-end tests for the regcorr() interface and the robustness fixes:
# convergence reporting, bootstrap filtering, and NA alignment.

make_normal_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  x <- runif(n)
  rho <- plogis(0.5 + 0.2 * x)
  z1 <- rnorm(n)
  y1 <- z1
  y2 <- z1 * rho + rnorm(n) * sqrt(1 - rho^2)
  data.frame(y1 = y1, y2 = y2, x = x)
}

test_that("regcorr fits bivariate normal responses with working methods", {
  dat <- make_normal_data()
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 30)

  expect_true(isTRUE(fit$converged))
  expect_equal(names(coef(fit)), c("(Intercept)", "x"))
  expect_length(fitted(fit), nrow(dat))
  expect_equal(nobs(fit), nrow(dat))
  expect_equal(dim(vcov(fit)), c(2, 2))
  expect_true(fit$nboot.valid > 20)

  pr <- predict(fit, newdata = data.frame(x = c(0.1, 0.5, 0.9)))
  expect_length(pr, 3)
  expect_true(all(pr > 0 & pr < 1))

  s <- summary(fit)
  expect_true(isTRUE(s$converged))
  expect_equal(s$nboot.valid, fit$nboot.valid)
})

test_that("predict returns NA for rows with missing covariates", {
  dat <- make_normal_data()
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 20)

  pr <- predict(fit, newdata = data.frame(x = c(0.2, NA, 0.8)))
  expect_length(pr, 3)
  expect_true(is.na(pr[2]))
  expect_false(anyNA(pr[c(1, 3)]))
})

test_that("predict handles NA factor levels", {
  set.seed(5)
  n <- 120
  g <- factor(sample(c("a", "b"), n, replace = TRUE))
  x <- runif(n)
  rho <- plogis(0.3 + 0.5 * (g == "b") + 0.2 * x)
  z1 <- rnorm(n)
  y1 <- z1
  y2 <- z1 * rho + rnorm(n) * sqrt(1 - rho^2)
  dat <- data.frame(y1 = y1, y2 = y2, g = g, x = x)

  # some bootstrap resamples may be unidentifiable (all one factor level)
  fit <- suppressWarnings(regcorr(cbind(y1, y2) ~ g + x, data = dat, nboot = 20))
  pr <- predict(fit, newdata = data.frame(g = factor(c("a", NA, "b")),
                                          x = c(0.3, 0.3, 0.7)))
  expect_length(pr, 3)
  expect_true(is.na(pr[2]))
})

test_that("missing values in the training data are dropped and recorded", {
  dat <- make_normal_data()
  dat$x[10] <- NA

  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 20)
  expect_equal(nobs(fit), nrow(dat) - 1)
  expect_false(is.null(fit$na.action))
  expect_length(fit$na.action, 1)

  # keeping NAs is not supported: clear error instead of silent garbage
  expect_error(
    regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 0,
            na.action = stats::na.pass),
    "Missing values"
  )
})

test_that("nboot = 0 yields NA vcov and no bootstrap", {
  dat <- make_normal_data(n = 100, seed = 7)
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 0)

  expect_true(all(is.na(vcov(fit))))
  expect_equal(fit$nboot.valid, 0)
})

test_that("bootstrap discards non-converged replications in degenerate binary data", {
  set.seed(11)
  n <- 40
  x2 <- runif(n)
  p1 <- plogis(1.2 * x2)
  p2 <- plogis(1.2 * x2)
  y1 <- rbinom(n, 1, p1)
  y2 <- rbinom(n, 1, p2)

  # the fit itself may also warn about non-convergence; collect all warnings
  warnings_seen <- character()
  fitb <- withCallingHandlers(
    regcorr(cbind(y1, y2) ~ x2, nboot = 60),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("discarded", warnings_seen)))
  expect_true(fitb$nboot.valid < 60)
  expect_true(fitb$nboot.valid > 2)   # enough replications to estimate vcov
  expect_true(is.matrix(vcov(fitb)))
  expect_equal(dim(vcov(fitb)), c(2, 2))
})

test_that("NR fit reports convergence status", {
  dat <- make_normal_data(n = 100, seed = 1)
  fit <- regcorr(cbind(y1, y2) ~ x, data = dat, nboot = 0)
  expect_true(isTRUE(fit$converged))
  # print output flags non-convergence when it occurs
  expect_output(print(fit), "Iterations")
})
