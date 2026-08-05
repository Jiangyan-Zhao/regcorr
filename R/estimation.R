#' Estimate beta for Bivariate Normal responses using Newton Raphson
#'
#' @param Y n by 2 matrix, paired responses.
#' @param X n by p matrix, covariate matrix including first column of ones.
#' @param betaIni Initial estimate of beta.
#' @param link Indicator of link function ("1" = logistic, "2" = tanh).
#' @return A list containing betaCurrent, numIter, restart, and converged
#'   (whether the iterations converged within the iteration/restart limits).
#' @importFrom stats lm.fit
NRfitBivNormal <- function(Y, X, betaIni, link)
{
  # return new estimate of beta using Newton Raphson method
  #
  # input:
  #        X:   n by p, covariate matrix including first column of ones
  #        Y:   n by 2, paired responses
  #    beta0:   initial estimate of beta


  # initial parameters
  numSample <- nrow(Y)
  p <- ncol(X)

  TOL <- 0.01
  maxIter <- 10
  maxRestart <- 10

  residual <- lm.fit(X, Y)$residuals                  # residuals of linear reg. on the full design matrix
  sigmaHat <- sqrt(colMeans(residual^2))              # std of y1 and y2
  Ytilde <- residual %*% diag(1 / sigmaHat)           # (y-muHat)/sigma
  T1 <- rowSums(Ytilde^2)                             # n by 1
  T2 <- apply(Ytilde, 1, prod)                        # n by 1

  betaCurrent <- betaIni
  numIter <- 0
  restart <- 0
  converged <- FALSE

  # main loop for Newton-Raphson
  while (numIter < maxIter && restart < maxRestart) {
    numIter <- numIter + 1
    betaPrev <- betaCurrent

    switch(link,
           "1" = { # logistic link
             rho <- logistic(X %*% betaCurrent)                    # n by 1, corr. coef.
             a <- -rho^2 / (1 + rho)^2 / (1 - rho)
             b <- rho * (1 + rho^2) / (1 + rho)^2 / (1 - rho)
             c <- rho^2 / (1 + rho)
             u <- -rho^2 * (2 - rho + rho^2) / (1 + rho)^3 / (1 - rho)
             v <- -rho * (-1 + rho - 5 * rho^2 + rho^3) / (1 + rho)^3 / (1 - rho)
             w <- rho^2 * (2 + rho) * (1 - rho) / (1 + rho)^2
           },
           "2" = { # tanh link
             rho <- tanh(X %*% betaCurrent)                        # n by 1, corr. coef.
             a <- -rho / (1 - rho) / (1 + rho)
             b <- (1 + rho^2) / (1 - rho) / (1 + rho)
             c <- rho
             u <- -(1 + rho^2) / (1 - rho) / (1 + rho)
             v <- 4 * rho / (1 - rho) / (1 + rho)
             w <- (1 - rho) * (1 + rho)
           }
    ) # end of switch

    abc <- a * T1 + b * T2 + c                         # n by 1
    uvw <- u * T1 + v * T2 + w                         # n by 1

    score <- t(X) %*% abc                              # p by 1
    H <- crossprod(X, X * as.vector(uvw))               # p by p

    ill_posed <- any(!is.finite(betaCurrent)) || any(!is.finite(score)) ||
      any(!is.finite(H)) || sqrt(sum(betaCurrent^2)) > 10 || kappa(H) > 10000
    if (ill_posed) {                                   # ill-posed Hessian/score
      # restart from a jittered starting value
      restart <- restart + 1
      numIter <- 0
      betaCurrent <- betaIni + stats::rnorm(p, 0, 0.1)
    } else {                                           # regular matrix
      step <- solve(H, score)
      betaCurrent <- betaCurrent - step
      if (sum((betaCurrent - betaPrev)^2) <= TOL) {    # converged
        converged <- TRUE
        break
      }
    }
  } # end of while

  list(betaCurrent = as.vector(betaCurrent), numIter = numIter,
       restart = restart, converged = converged)
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
  if (link == "1") {                    # logistic link
    hp  <- rho * (1 - rho)              # h'(eta) = rho*(1-rho)
    hpp <- rho * (1 - rho) * (1 - 2*rho)  # h''(eta)
  } else {                              # tanh link
    hp  <- 1 - rho^2                    # h'(eta) = 1 - rho^2
    hpp <- -2 * rho * (1 - rho^2)       # h''(eta)
  }
  denom <- c + d * rho
  list(
    score_weight   = d * hp / denom,
    hessian_weight = -d^2 * hp^2 / denom^2 + d * hpp / denom
  )
}

#' Estimate beta for Bivariate Bernoulli responses using Newton Raphson
#'
#' @param Y n by 2 matrix, paired responses.
#' @param X n by p matrix, covariate matrix including first column of ones.
#' @param beta0 Initial estimate of beta.
#' @param link Indicator of link function ("1" = logistic, "2" = tanh).
#' @param warn Whether to warn when the marginal logistic fits hit fitted
#'   probabilities at (or near) 0/1, indicating perfect or quasi-perfect
#'   separation. Suppressed during bootstrap replications.
#' @return A list containing betaCurrent, numIter, restart, and converged
#'   (whether the iterations converged within the iteration/restart limits).
#' @importFrom stats glm.fit binomial
NRfitBivBernoulli <- function(Y, X, beta0, link, warn = FALSE)
{
  # return new estimate of beta using Newton Raphson method
  #
  # input:
  #        X:   n by p, covariate matrix including first column of ones
  #        Y:   n by 2, paired responses
  #    beta0:   initial estimate of beta
  #     link:   indicator of link function, "1" = logistic, "2" = tanh

  # initial parameters
  p <- ncol(X)
  n <- nrow(Y)
  TOL <- 0.01
  maxIter <- 10
  maxRestart <- 10

  # get marginal pHat via logistic regressions on the full design matrix
  eta1Hat <- suppressWarnings(
    glm.fit(X, Y[, 1], family = binomial(link = "logit"))$coefficients
  )
  eta2Hat <- suppressWarnings(
    glm.fit(X, Y[, 2], family = binomial(link = "logit"))$coefficients
  )

  p1Hat <- as.matrix(logistic(X %*% eta1Hat), ncol = 1)  # marginal est of success prob for y1
  p2Hat <- as.matrix(logistic(X %*% eta2Hat), ncol = 1)  # marginal est of success prob for y2

  if (warn && (any(p1Hat <= 1e-8) || any(p1Hat >= 1 - 1e-8) ||
               any(p2Hat <= 1e-8) || any(p2Hat >= 1 - 1e-8))) {
    warning("Marginal logistic fits produced fitted probabilities near 0 or 1 ",
            "(perfect or quasi-perfect separation); the correlation parameter ",
            "may be weakly identified.",
            call. = FALSE)
  }

  I00 <- (1 - Y[, 1]) * (1 - Y[, 2])   # n by 1
  I01 <- (1 - Y[, 1]) * Y[, 2]         # n by 1
  I10 <- Y[, 1] * (1 - Y[, 2])         # n by 1
  I11 <- Y[, 1] * Y[, 2]               # n by 1

  c11 <- p1Hat * p2Hat
  c10 <- p1Hat - p1Hat * p2Hat
  c01 <- p2Hat - p1Hat * p2Hat
  c00 <- 1 - p1Hat - p2Hat + p1Hat * p2Hat

  d00 <- d11 <- sqrt(c11 * c00)
  d10 <- d01 <- -d11

  betaCurrent <- beta0
  numIter <- 0
  restart <- 0
  converged <- FALSE

  while (numIter < maxIter && restart < maxRestart) { # when beta not converge
    numIter <- numIter + 1
    betaPrev <- betaCurrent
    switch(link,
           "1" = { # logistic link
             rho <- logistic(X %*% betaCurrent)     # n by 1, corr. coef.
           },
           "2" = { # tanh link
             rho <- tanh(X %*% betaCurrent)
           }
    ) # end of switch

    # score and Hessian weights for each cell of log(c_ab + d_ab*rho)
    a00 <- .bernoulli_derivatives(c00, d00, rho, link)
    a01 <- .bernoulli_derivatives(c01, d01, rho, link)
    a10 <- .bernoulli_derivatives(c10, d10, rho, link)
    a11 <- .bernoulli_derivatives(c11, d11, rho, link)
    e00 <- a00$score_weight;  f00 <- a00$hessian_weight
    e01 <- a01$score_weight;  f01 <- a01$hessian_weight
    e10 <- a10$score_weight;  f10 <- a10$hessian_weight
    e11 <- a11$score_weight;  f11 <- a11$hessian_weight

    Z <- I00 * e00 + I01 * e01 + I10 * e10 + I11 * e11   # n by 1
    score <- t(X) %*% Z                                  # p by 1

    w <- I00 * f00 + I01 * f01 + I10 * f10 + I11 * f11   # n by 1

    H <- crossprod(X, X * as.vector(w))                   # p by p

    ill_posed <- any(!is.finite(betaCurrent)) || any(!is.finite(score)) ||
      any(!is.finite(H)) || sqrt(sum(betaCurrent^2)) > 10 || kappa(H) > 10000
    if (ill_posed) {                                     # ill-posed Hessian/score
      # restart from a jittered starting value
      restart <- restart + 1
      numIter <- 0
      betaCurrent <- beta0 + stats::rnorm(p, 0, 0.1)
    } else {                                             # regular matrix
      step <- solve(H, score)
      betaCurrent <- betaCurrent - step
      if (sum((betaCurrent - betaPrev)^2) <= TOL) {      # converged
        converged <- TRUE
        break
      }
    }
  } # end of while

  list(betaCurrent = as.vector(betaCurrent), numIter = numIter,
       restart = restart, converged = converged)
}
