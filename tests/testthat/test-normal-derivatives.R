# Numerical validation of the normal score and Hessian after the marginal
# estimates have been fixed, for both correlation links.

test_that("normal observation-level derivatives match finite differences", {
  set.seed(2718)
  eta <- seq(-1.2, 1.2, length.out = 11)
  T1 <- runif(length(eta), 0.5, 4)
  T2 <- runif(length(eta), -0.6, 0.6)

  for (link in c("1", "2")) {
    rho <- .correlation_link(eta, link)
    analytic <- .normal_derivatives(T1, T2, rho, link)
    loglik_eta <- function(value) {
      .normal_loglik_contributions(
        T1, T2, .correlation_link(value, link)
      )
    }

    score_numerical <- central_derivative(loglik_eta, eta)
    hessian_numerical <- central_second_derivative(loglik_eta, eta)

    expect_equal(
      analytic$score_weight, score_numerical,
      tolerance = 1e-6
    )
    expect_equal(
      analytic$hessian_weight, hessian_numerical,
      tolerance = 2e-4
    )
  }
})

test_that("assembled normal score and Hessian match finite differences", {
  set.seed(31415)
  n <- 90
  X <- cbind(1, x1 = runif(n, -1, 1), x2 = rnorm(n, sd = 0.5))
  rho_data <- tanh(X %*% c(0.1, 0.25, -0.15))
  z1 <- rnorm(n)
  Y <- cbind(
    y1 = 0.3 + X[, "x1"] + z1,
    y2 = -0.2 + 0.5 * X[, "x2"] +
      rho_data * z1 + sqrt(1 - rho_data^2) * rnorm(n)
  )
  marginal <- .normal_marginal_statistics(Y, X)
  expect_true(marginal$ok)

  points <- list(
    `1` = c(-0.4, 0.2, -0.1),
    `2` = c(0.1, 0.2, -0.1)
  )
  for (link in c("1", "2")) {
    beta <- points[[link]]
    analytic <- .normal_objective(
      beta, X, marginal$T1, marginal$T2, link
    )
    expect_true(analytic$valid)

    loglik <- function(value) {
      .normal_objective(
        value, X, marginal$T1, marginal$T2, link
      )$loglik
    }
    score_numerical <- central_gradient(loglik, beta)
    hessian_numerical <- central_hessian(loglik, beta)

    expect_equal(analytic$score, score_numerical, tolerance = 1e-5)
    expect_equal(
      unname(analytic$hessian), unname(hessian_numerical),
      tolerance = 1e-3
    )
  }
})
