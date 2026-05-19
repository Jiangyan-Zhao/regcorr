rbinary <-function(n,p,rho)
{
# generate n bivariate binary row vectors with mean p and correlation rho
# input:
#     n: number of rows
#     p: 1 by 2 mean vector of bivariate variables
#     rho: correlation of bivariate variables
# reference:
#   Qaqish, B. F. (2003). A family of multivariate binary distributions for
#      simulating correlated binary variables with specified marginal means
#      and correlations. Biometrika 90, 455-463.
#
# example:
# n=10; p=c(0.4,0.5); rho=0.5
# y <- rbinary(10,c(0.4,0.5),0.5)

b=rho*sqrt(p[2]*(1-p[2])/p[1]/(1-p[1]))
Y=matrix(0,n,2)
Y[,1]=rbinom(n,1,p[1])
Y[,2]=rbinom(n,1,p[2]+b*(Y[,1]-p[1]))
return(Y)
}
