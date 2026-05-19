#' Subroutine to test the significance of individual parameters
#'
#' @param numSample Sample size.
#' @param p Number of covariates.
#' @param link Link function.
#' @param model Model type (1: biv normal; 2: biv Bernoulli).
#' @param betaTrue True beta.
#' @param betaIni Initial beta.
#' @param eta1True True eta1.
#' @param eta2True True eta2.
#' @param numSimu Number of simulations.
#' @param numBoot Number of bootstrap iterations.
#' @return A list containing RMSE, ConsistRate, and testing power.
#' @importFrom stats pchisq pnorm cov
#' @export
subRoutineTest<-function(numSample,p,link,model,betaTrue,betaIni,eta1True,eta2True,numSimu,numBoot)
{ # return powers of testing for significance of individual parameters

  pval=matrix(0,numSimu,p+2)       # globe and indi. pvals
  cBeta=matrix(0,numSimu,p+1)      # storage for betaHat
  cConsistRate=matrix(0,numSimu,1) # storage for consistency rate

  for(iSimu in 1:numSimu){

    ## (1) estimation

    # generate data and estimate beta
    switch(model, # 1: biv normal; 2: biv Bernoulli
           "1" = {
             try_res <- try(NRfitBivNormal(YBoot,XBoot,betaIni,link)$betaCurrent, silent = TRUE)
             if(!inherits(try_res, "try-error")) betaHatBoot[iBoot,] = try_res else betaHatBoot[iBoot,] = betaIni
           },
           "2" = {
             try_res <- try(NRfitBivBernoulli(YBoot,XBoot,betaIni,link)$betaCurrent, silent = TRUE)
             if(!inherits(try_res, "try-error")) betaHatBoot[iBoot,] = try_res else betaHatBoot[iBoot,] = betaIni
           }
    ) # end of switch
    cBeta[iSimu,]=betaHat    # collect estimates

    # estimate consistency rate
    switch(link, # 1: logistic; 2??tanh
           "1" = {rhoHat=as.matrix(logistic(data$X%*%cBeta[iSimu,]),ncol=1)},
           "2" = {rhoHat=as.matrix(tanh(data$X%*%cBeta[iSimu,]),ncol=1)}
    ) # end of switch
    cConsistRate[iSimu]=mean(rhoHat*data$rho>0)

    ## (2) testing

    # Bootstrap estimation for cov of betaHat

    betaHatBoot=matrix(0,numBoot,p+1)
    for(iBoot in 1:numBoot){
      indexBoot=sample(numSample,replace=TRUE)
      YBoot=data$Y[indexBoot,];  XBoot=data$X[indexBoot,] # boot data

      switch(model, # 1: biv normal; 2: biv Bernoulli
             "1" = {betaHatBoot[iBoot,]=NRfitBivNormal(YBoot,XBoot,betaIni,link)$betaCurrent   },
             "2" = {betaHatBoot[iBoot,]=NRfitBivBernoulli(YBoot,XBoot,betaIni,link)$betaCurrent}
      ) # end of switch
    } # end of bootstrap

    badBoot=which(rowSums(abs(betaHatBoot-matrix(betaIni,1)[rep(1,numBoot),]))<0.01)
    if (length(badBoot)>0){              # clear bad boot if any
      betaHatBoot=betaHatBoot[-badBoot,] # remove bad bootstrap items
      print(paste("# of bad boots is ", length(badBoot)))
    } # end of if

    # compute p-values
    if (length(badBoot)<450) { # acceptable data for bootstrap
      covBetaHat=cov(betaHatBoot)
      pval[iSimu,1]=1-pchisq(t(betaHat[-1])%*%
                               solve(covBetaHat[-1,-1])%*%betaHat[-1],p) # global p-value
      pval[iSimu,2:(p+2)]=2*(1-pnorm(abs(betaHat)/sqrt(diag(covBetaHat)))) # indi. p-value
    } else                 # bad data for bootstrap
    {iSimu=max(iSimu-1,1)} # re-do this iSimu

  } # end of simu

  # report
  RMSE=sqrt(mean(rowSums((cBeta-matrix(betaTrue,1)[rep(1,dim(cBeta)[1]),])^2)))
  ConsistRate=mean(cConsistRate)
  power=colMeans(pval<0.05,na.rm=TRUE)

  return(list(RMSE=RMSE,ConsistRate=ConsistRate,power=power))
}
