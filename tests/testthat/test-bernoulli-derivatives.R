# Numerical-derivative validation of the Bernoulli score and Hessian weights.
#
# The log-likelihood of the bivariate Bernoulli model is
#   ell(beta) = sum_i sum_ab I_ab * log(c_ab + d_ab * rho_i),
# with rho_i = h(X_i beta).  Its score and Hessian are
#   score   = t(X) %*% Z,               Z  = sum_ab I_ab * e_ab,
#   Hessian = t(X) %*% diag(w) %*% X,   w  = sum_ab I_ab * f_ab,
# where
#   e_ab = d_ab * h'(eta) / (c_ab + d_ab * rho)
#   f_ab = -d_ab^2 * (h'(eta))^2 / (c_ab + d_ab * rho)^2
#          + d_ab * h''(eta) / (c_ab + d_ab * rho)
# with
#   h'   = rho*(1-rho),            h'' = rho*(1-rho)*(1-2*rho)   (logistic)
#   h'   = 1 - rho^2,              h'' = -2*rho*(1-rho^2)         (tanh)
#
# These analytic expressions are checked against central finite
# differences of ell(beta).

cell_loglik <- function(eta, c, d, link) {
  rho <- if (link == "1") logistic(eta) else tanh(eta)
  log(c + d * rho)
}

num_deriv1 <- function(f, x, h = 1e-6) {
  (f(x + h) - f(x - h)) / (2 * h)
}

num_deriv2 <- function(f, x, h = 1e-4) {
  (f(x + h) - 2 * f(x) + f(x - h)) / h^2
}

num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x
    xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(1))
}

num_hess <- function(f, x, h = 1e-4) {
  p <- length(x)
  H <- matrix(0, p, p)
  for (j in seq_len(p)) {
    xp <- x
    xm <- x
    xp[j] <- x[j] + h
    xm[j] <- x[j] - h
    H[j, j] <- (f(xp) - 2 * f(x) + f(xm)) / h^2
    if (j < p) {
      for (k in (j + 1):p) {
        xpp <- x; xpm <- x; xmp <- x; xmm <- x
        xpp[j] <- x[j] + h; xpp[k] <- x[k] + h
        xpm[j] <- x[j] + h; xpm[k] <- x[k] - h
        xmp[j] <- x[j] - h; xmp[k] <- x[k] + h
        xmm[j] <- x[j] - h; xmm[k] <- x[k] - h
        H[j, k] <- H[k, j] <- (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h^2)
      }
    }
  }
  H
}

test_that("single-cell score and Hessian weights match numerical derivatives", {
  set.seed(7)
  for (link in c("1", "2")) {
    eta <- seq(-2.5, 2.5, length.out = 9)
    # keep c + d*rho > 0 for all cells so log() is well defined
    c <- runif(9, 0.3, 0.7)
    d <- runif(9, -0.3, 0.3)
    rho <- if (link == "1") logistic(eta) else tanh(eta)

    deriv <- .bernoulli_derivatives(c, d, rho, link)

    e_num <- num_deriv1(function(x) cell_loglik(x, c, d, link), eta)
    f_num <- num_deriv2(function(x) cell_loglik(x, c, d, link), eta)

    expect_equal(deriv$score_weight, e_num, tolerance = 1e-6)
    expect_equal(deriv$hessian_weight, f_num, tolerance = 1e-4)
  }
})

test_that("assembled Bernoulli score and Hessian match numerical derivatives of the log-likelihood", {
  set.seed(99)
  for (link in c("1", "2")) {
    beta_true <- c(0.2, 0.4, -0.3)
    dat <- genDataBB(
      numSample = 100, p = 2,
      betaTrue = beta_true,
      eta1True = c(0.1, 0.5, 0.2),
      eta2True = c(-0.2, 0.3, 0.5),
      link = link
    )
    X <- dat$X
    Y <- dat$Y

    eta1Hat <- glm.fit(X, Y[, 1], family = binomial("logit"))$coefficients
    eta2Hat <- glm.fit(X, Y[, 2], family = binomial("logit"))$coefficients
    p1 <- logistic(X %*% eta1Hat)
    p2 <- logistic(X %*% eta2Hat)

    c11 <- p1 * p2
    c10 <- p1 - p1 * p2
    c01 <- p2 - p1 * p2
    c00 <- 1 - p1 - p2 + p1 * p2
    d00 <- d11 <- sqrt(c11 * c00)
    d10 <- d01 <- -d11

    I00 <- (1 - Y[, 1]) * (1 - Y[, 2])
    I01 <- (1 - Y[, 1]) * Y[, 2]
    I10 <- Y[, 1] * (1 - Y[, 2])
    I11 <- Y[, 1] * Y[, 2]

    rho <- as.vector(if (link == "1") logistic(X %*% beta_true) else tanh(X %*% beta_true))

    a00 <- .bernoulli_derivatives(c00, d00, rho, link)
    a01 <- .bernoulli_derivatives(c01, d01, rho, link)
    a10 <- .bernoulli_derivatives(c10, d10, rho, link)
    a11 <- .bernoulli_derivatives(c11, d11, rho, link)

    Z <- I00 * a00$score_weight + I01 * a01$score_weight +
      I10 * a10$score_weight + I11 * a11$score_weight
    w <- I00 * a00$hessian_weight + I01 * a01$hessian_weight +
      I10 * a10$hessian_weight + I11 * a11$hessian_weight

    g_analytic <- as.vector(t(X) %*% Z)
    H_analytic <- t(X) %*% (X * as.vector(w))

    loglik <- function(b) {
      eta <- X %*% b
      rho_ll <- as.vector(if (link == "1") logistic(eta) else tanh(eta))
      sum(I00 * log(c00 + d00 * rho_ll) + I01 * log(c01 + d01 * rho_ll) +
            I10 * log(c10 + d10 * rho_ll) + I11 * log(c11 + d11 * rho_ll))
    }

    g_num <- num_grad(loglik, beta_true)
    H_num <- num_hess(loglik, beta_true)

    expect_equal(g_analytic, g_num, tolerance = 1e-5)
    expect_equal(H_analytic, H_num, tolerance = 1e-3)
  }
})
