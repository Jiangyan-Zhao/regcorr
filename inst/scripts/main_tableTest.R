library(regcorr)

numSimu=1000; numBoot=500;
cSample=c(50,100);lenCS=length(cSample)

tableTest=NULL #

system.time({

for (iModel in 2:2){   # 1: biv normal; 2:biv Bernoulli
  for (iLink in 1:1) { # 1: logistic; 2:tanh

    res_model_link=NULL # clear storage

    for(iCova in 2:3){ # number of covariates

    print(c(iModel,iLink,iCova))

      # assign true parameters
      eta1True=eta2True=as.vector(rep(0,iCova+1))   #  for marginal y1, y2
      switch(iLink, # 1: logistic; 2:tanh
             "1" = {betaTrue=as.vector(c(0.25,1,rep(0,iCova-1)))  # by latex set up
             betaIni=as.vector(c(0.25,rep(0,iCova)))      },
             "2" = {betaTrue=as.vector(c(0.25,-1,rep(0,iCova-1))) # by latex set up
             betaIni=as.vector(c(0.25,rep(0,iCova)))     }
      ) # end of switch link

      res_model_link_p=NULL        # clear storage
      for (iNumSample in cSample){ # results under fixed model, link and p
          res=subRoutineTest(iNumSample,iCova,iLink,iModel,
                             betaTrue,betaIni,eta1True,eta2True,numSimu,numBoot)
          res_model_link_p=rbind(res_model_link_p, # collect over numSample
                                 c(res$RMSE,         # RMSE
                                 res$ConsistRate,  # consistency rate
                                 res$power[-2]))    # pvals except beta-0
      } # end of iNumSample

      res_model_link=cbind(res_model_link, res_model_link_p)  # collect over cova

    } # end of iCova

    tableTest=rbind(tableTest,cbind(rep(iModel,lenCS),
                                    rep(iLink,lenCS),
                                    cSample,res_model_link)) # collect over model and link
  } # end of iLink
}   # end of iModel
})  # end of system.time

tableTest=as.data.frame(tableTest)
names(tableTest)=c("model","link","n","RMSE","CR","bG","b1","b2",
                                      "RMSE","CR","bG","b1","b2","b3")
tableTest
if (!dir.exists("inst/extdata")) dir.create("inst/extdata", recursive = TRUE)
write.csv(tableTest, file = "inst/extdata/tableTestPower.csv", row.names = FALSE)

