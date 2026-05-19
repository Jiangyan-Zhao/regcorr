NRfitBivBernoulli <- function(Y,X,beta0,link)
{
  # return new estimate of beta using Newton Raphson method
  #
  # input:
  #        X:   n by p, covariate matrix including first column of ones
  #        Y:   n by 2, paired responses
  #    beta0:   initial estimate of beta
  #     link:   indicator of link function, "1" = logistic, "2" = tanh
  
  # initial parameters
  p = ncol(X)
  n = nrow(Y)
  TOL=0.01
  
  # get marginal pHat
  switch(p-1, # num of covariates
         "1" = {  eta1Hat <- as.matrix(glm(Y[,1] ~ X[,2],family = binomial(link="logit"))$coefficients,ncol=1)
                  eta2Hat <- as.matrix(glm(Y[,2] ~ X[,2],family = binomial(link="logit"))$coefficients,ncol=1)
         },
         "2" = {  eta1Hat <- as.matrix(glm(Y[,1] ~ X[,2]+X[,3],family = binomial(link="logit"))$coefficients,ncol=1)
                  eta2Hat <- as.matrix(glm(Y[,2] ~ X[,2]+X[,3],family = binomial(link="logit"))$coefficients,ncol=1)
         },
         "3" = {  eta1Hat <- as.matrix(glm(Y[,1] ~ X[,2]+X[,3]+X[,4],family = binomial(link="logit"))$coefficients,ncol=1)
                  eta2Hat <- as.matrix(glm(Y[,2] ~ X[,2]+X[,3]+X[,4],family = binomial(link="logit"))$coefficients,ncol=1)
         }
  ) # end of switch
  
  p1Hat=as.matrix(logistic(X%*%eta1Hat),ncol=1)  # marginal est of success prob for y1, numSample by 1
  p2Hat=as.matrix(logistic(X%*%eta2Hat),ncol=1)  # marginal est of success prob for y2, numSample by 1
  
  # 
  I00=(1-Y[,1])*(1-Y[,2])  # n by 1
  I01=(1-Y[,1])*Y[,2]      # n by 1
  I10=Y[,1]*(1-Y[,2])      # n by 1
  I11=Y[,1]*Y[,2]          # n by 1
  
  c11=p1Hat*p2Hat
  c10=p1Hat-p1Hat*p2Hat
  c01=p2Hat-p1Hat*p2Hat
  c00=1-p1Hat-p2Hat+p1Hat*p2Hat
  
  d00=d11=sqrt(c11*c00)
  d10=d01=-d11
  
  # initial beta
  betaPrev=runif(p)     # dummy start
  betaCurrent=beta0
  restart=numIter=0
 
  
  while(sum((betaCurrent-betaPrev)^2)>TOL & numIter<10 & restart<10) # when beta not converge
  { 
    numIter=numIter+1
    betaPrev=betaCurrent
    switch(link,
           "1" = { # logistic link
            rho=logistic(X%*%betaCurrent)     # n by 1, corr. coef.
            e00=d00*rho*(1-rho)/(c00+d00*rho) # n by 1
            e01=d01*rho*(1-rho)/(c01+d01*rho) # n by 1
            e10=d10*rho*(1-rho)/(c10+d10*rho) # n by 1
            e11=d11*rho*(1-rho)/(c11+d11*rho) # n by 1
            f00=(-d00^2+d00*(c00+d00*rho)*(1-2*rho))*rho^2*(1-rho)^2/(c00+d00*rho)^2 # n by 1
            f01=(-d01^2+d01*(c01+d01*rho)*(1-2*rho))*rho^2*(1-rho)^2/(c01+d01*rho)^2 # n by 1
            f10=(-d10^2+d10*(c10+d10*rho)*(1-2*rho))*rho^2*(1-rho)^2/(c10+d10*rho)^2 # n by 1
            f11=(-d11^2+d11*(c11+d11*rho)*(1-2*rho))*rho^2*(1-rho)^2/(c11+d11*rho)^2 # n by 1
           },
           "2" = { # tanh link 
             rho=tanh(X%*%betaCurrent)
             e00=d00*(1+rho)*(1-rho)/(c00+d00*rho) # n by 1
             e01=d01*(1+rho)*(1-rho)/(c01+d01*rho) # n by 1
             e10=d10*(1+rho)*(1-rho)/(c10+d10*rho) # n by 1
             e11=d11*(1+rho)*(1-rho)/(c11+d11*rho) # n by 1
             f00=(-d00^2-d00*(c00+d00*rho)*2*rho)*(1+rho)^2*(1-rho)^2/(c00+d00*rho)^2 # n by 1
             f01=(-d01^2-d01*(c01+d01*rho)*2*rho)*(1+rho)^2*(1-rho)^2/(c01+d01*rho)^2 # n by 1
             f10=(-d10^2-d10*(c10+d10*rho)*2*rho)*(1+rho)^2*(1-rho)^2/(c10+d10*rho)^2 # n by 1
             f11=(-d11^2-d11*(c11+d11*rho)*2*rho)*(1+rho)^2*(1-rho)^2/(c11+d11*rho)^2 # n by 1
             }
    ) # end of switch
    
    Z=I00*e00+I01*e01+I10*e10+I11*e11                     # n by 1
    score=t(X) %*% Z                                      # p by 1
    
    w=I00*f00+I01*f01+I10*f10+I11*f11                     # n by 1
    
    H=t(X) %*% diag(as.vector(w)) %*% X                   # p by p

    if (sqrt(sum(betaCurrent^2))>10 || kappa(H)>10000 ) { # ill posed matrix
      betaPrev=runif(p); betaCurrent=beta0; numIter=1
      restart=restart+1 
     } else {                                          # regular matrix
     invH=solve(H)                                     # p by p, inverse Hessian
     betaCurrent <- betaCurrent-invH%*%score           # update beta
    }  # end of if 

  } # end of while
  #if (restart>1) {print(paste("restart=",restart))}
  
return(result=list(betaCurrent=betaCurrent,numIter=numIter,restart=restart))
}
