# genDataBB must keep exactly the rows whose correlation is feasible for the
# given marginal probabilities, i.e. the row-wise bounds of Dufera et al.
# (2023), Section 6 (Qaqish, 2003):
#   max(-psi1*psi2, -(psi1*psi2)^-1) <= rho <= min(psi1/psi2, psi2/psi1)
# with psi_j = sqrt(p_j / (1 - p_j)). This guards against the previous
# column-wise (apply MARGIN = 2) implementation that dropped valid rows.

test_that("genDataBB keeps exactly the rows satisfying the Qaqish bounds", {
  for (link in c("1", "2")) {
    # replicate genDataBB's own X generation (the first RNG call after the
    # seed) so the feasibility check below is computed on the same rows
    set.seed(1)
    X_full <- cbind(rep(1, 500), matrix(runif(500 * 2), 500))
    eta1 <- c(0, 2, -1)
    eta2 <- c(0, -1.5, 1)
    beta <- c(0.25, 1, 0)

    rho_full <- as.vector(if (link == "1") logistic(X_full %*% beta)
                          else tanh(X_full %*% beta))
    p1 <- as.vector(logistic(X_full %*% eta1))
    p2 <- as.vector(logistic(X_full %*% eta2))
    psi1 <- sqrt(p1 / (1 - p1))
    psi2 <- sqrt(p2 / (1 - p2))
    feasible <- pmax(-psi1 * psi2, -1 / (psi1 * psi2)) <= rho_full &
      rho_full <= pmin(psi1 / psi2, psi2 / psi1)

    set.seed(1)
    dat <- genDataBB(500, 2, betaTrue = beta,
                     eta1True = eta1, eta2True = eta2, link = link)

    expect_equal(nrow(dat$X), sum(feasible))

    # all kept rows satisfy the bounds (and hence have positive cell
    # probabilities, so the correlation is well defined)
    p1_keep <- as.vector(logistic(dat$X %*% eta1))
    p2_keep <- as.vector(logistic(dat$X %*% eta2))
    psi1_keep <- sqrt(p1_keep / (1 - p1_keep))
    psi2_keep <- sqrt(p2_keep / (1 - p2_keep))
    lb <- pmax(-psi1_keep * psi2_keep, -1 / (psi1_keep * psi2_keep))
    ub <- pmin(psi1_keep / psi2_keep, psi2_keep / psi1_keep)
    expect_true(all(lb <= dat$rho & dat$rho <= ub))
  }
})

test_that("genDataBB produces valid bivariate binary data", {
  set.seed(2)
  for (link in c("1", "2")) {
    dat <- genDataBB(500, 2, betaTrue = c(0.25, 1, 0),
                     eta1True = c(0, 1, -0.5), eta2True = c(0, -1, 0.5),
                     link = link)
    expect_true(nrow(dat$Y) > 0)
    expect_true(all(dat$Y %in% c(0, 1)))
    # empirical correlation should be consistent with the modelled rho; use
    # an absolute tolerance because expect_equal()'s tolerance is relative
    # (15% of mean rho ~ 0.5 allows only ~0.075 slack, tighter than the
    # sampling noise from ~300 retained rows)
    expect_lt(abs(cor(dat$Y[, 1], dat$Y[, 2]) - mean(dat$rho)), 0.15)
  }
})
