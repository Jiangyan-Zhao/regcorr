#' Subroutine to Test the Significance of Individual Parameters
#'
#' Conduct simulation studies for testing regression coefficients in
#' correlation regression models.
#'
#' @param numSample Sample size.
#' @param p Number of covariates.
#' @param link Link function. "1" for logistic and "2" for tanh.
#' @param model Model type. "1" for bivariate normal and "2" for bivariate Bernoulli.
#' @param betaTrue True beta coefficients.
#' @param betaIni Initial beta coefficients.
#' @param eta1True True eta1 coefficients.
#' @param eta2True True eta2 coefficients.
#' @param numSimu Number of simulation replications.
#' @param numBoot Number of bootstrap iterations.
#'
#' @return A list containing:
#' \describe{
#'   \item{RMSE}{Root mean squared error of the estimated coefficients.}
#'   \item{ConsistRate}{Consistency rate of estimated correlations.}
#'   \item{power}{Estimated testing power.}
#' }
#'
#' @examples
#' set.seed(123)
#'
#' res <- subRoutineTest(
#'   numSample = 20,
#'   p = 1,
#'   link = "1",
#'   model = "1",
#'   betaTrue = c(0.2, 0.1),
#'   betaIni = c(0, 0),
#'   eta1True = c(0, 0),
#'   eta2True = c(0, 0),
#'   numSimu = 1,
#'   numBoot = 5
#' )
#'
#' names(res)
#'
#' @importFrom stats pchisq pnorm cov
#' @export
subRoutineTest <- function(numSample, p, link, model,
                           betaTrue, betaIni,
                           eta1True, eta2True,
                           numSimu, numBoot)
{
  # return powers of testing for significance of individual parameters

  pval = matrix(0, numSimu, p + 2)       # global and individual p-values
  cBeta = matrix(0, numSimu, p + 1)      # storage for betaHat
  cConsistRate = matrix(0, numSimu, 1)   # storage for consistency rate

  for(iSimu in 1:numSimu){

    ## (1) estimation

    switch(model,
           "1" = {
             data = genDataBN(
               numSample,
               p,
               betaTrue,
               eta1True,
               eta2True,
               link
             )

             betaHat =
               NRfitBivNormal(
                 data$Y,
                 data$X,
                 betaIni,
                 link
               )$betaCurrent
           },

           "2" = {
             data = genDataBB(
               numSample,
               p,
               betaTrue,
               eta1True,
               eta2True,
               link
             )

             betaHat =
               NRfitBivBernoulli(
                 data$Y,
                 data$X,
                 betaIni,
                 link
               )$betaCurrent
           }
    )

    cBeta[iSimu, ] = betaHat

    ## consistency rate

    switch(link,
           "1" = {
             rhoHat =
               as.matrix(
                 logistic(data$X %*% cBeta[iSimu, ]),
                 ncol = 1
               )
           },

           "2" = {
             rhoHat =
               as.matrix(
                 tanh(data$X %*% cBeta[iSimu, ]),
                 ncol = 1
               )
           }
    )

    cConsistRate[iSimu] = mean(rhoHat * data$rho > 0)

    ## (2) bootstrap

    betaHatBoot = matrix(0, numBoot, p + 1)

    for(iBoot in 1:numBoot){

      indexBoot = sample(numSample,
                         replace = TRUE)

      YBoot = data$Y[indexBoot, ]
      XBoot = data$X[indexBoot, ]

      switch(model,

             "1" = {
               betaHatBoot[iBoot, ] =
                 NRfitBivNormal(
                   YBoot,
                   XBoot,
                   betaIni,
                   link
                 )$betaCurrent
             },

             "2" = {
               betaHatBoot[iBoot, ] =
                 NRfitBivBernoulli(
                   YBoot,
                   XBoot,
                   betaIni,
                   link
                 )$betaCurrent
             }
      )
    }

    badBoot =
      which(
        rowSums(
          abs(
            betaHatBoot -
              matrix(betaIni, 1)[rep(1, numBoot), ]
          )
        ) < 0.01
      )

    if(length(badBoot) > 0){
      betaHatBoot = betaHatBoot[-badBoot, , drop = FALSE]
    }

    ## compute p-values

    if(nrow(betaHatBoot) > p + 1 &&
       length(badBoot) < 450){

      covBetaHat = cov(betaHatBoot)

      pval[iSimu, 1] =
        1 - pchisq(
          t(betaHat[-1]) %*%
            solve(covBetaHat[-1, -1]) %*%
            betaHat[-1],
          p
        )

      pval[iSimu, 2:(p + 2)] =
        2 * (
          1 -
            pnorm(
              abs(betaHat) /
                sqrt(diag(covBetaHat))
            )
        )
    }
  }

  ## report

  RMSE =
    sqrt(
      mean(
        rowSums(
          (
            cBeta -
              matrix(betaTrue, 1)[rep(1, nrow(cBeta)), ]
          )^2
        )
      )
    )

  ConsistRate = mean(cConsistRate)

  power = colMeans(
    pval < 0.05,
    na.rm = TRUE
  )

  return(
    list(
      RMSE = RMSE,
      ConsistRate = ConsistRate,
      power = power
    )
  )
}
