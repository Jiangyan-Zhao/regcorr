genDataBN<-function(numSample,p,betaTrue,eta1True,eta2True,link)
{ # return data generated from biv normal
  
  # generate X, numSample by numCova+1, all uniform
  X=cbind(rep(1,numSample),matrix(runif(numSample*p),numSample)) 

    switch(link, # 1: logistic; 2:tanh
         "1" = {rho=as.matrix(logistic(X%*%betaTrue),ncol=1) },
         "2" = {rho=as.matrix(tanh(X%*%betaTrue),ncol=1) }
  ) # end of switch
  
  sigma1=1; sigma2=1
  # generate Y with eta1 and eta2 both zero vector
  Y=matrix(0,numSample,2)
  z1=rnorm(numSample)
  Y[,1]=sigma1*z1                                       # Y[,1]~N(0,sigma1^2)
  Y[,2]=sigma2*(z1*rho+rnorm(numSample)*sqrt(1-rho^2))  # Y[,2]~N(0,sigma2^2)

  return(list(X=X,Y=Y,rho=rho))
} 
