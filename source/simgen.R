
simple.setup = function(N,p,alpha,beta1,beta0,
                        st.mu=1,sy.mu=1)
{
  X = matrix(rnorm(N*p),N,p)
  ps = 1/(1+exp(-drop(X %*% alpha[-1] + alpha[1])))
  trt = rbinom(N,1,ps)
  
  pY1 =  1/(1+exp(-drop(X %*% beta1[-1]+beta1[1])))
  pY0 = 1/(1+exp(- drop(X %*% beta0[-1]+beta0[1])))
  pY = trt*pY1 + (1-trt)*pY0
  Y = rbinom(N,1,pY)

  S = cbind( trt*st.mu + rnorm(N),
                  Y*sy.mu + rnorm(N))
  
  return(list(Y = Y, trt = trt, 
                    X = X, S = S,
              ate = mean(pY1-pY0)))
}

sim.mix = function(N,p,alpha,beta1,beta0,
                       st.mu=1,sy.mu=1)
{
  X = matrix(rnorm(N*p),N,p)
  ps = 1/(1+exp(-drop(X %*% alpha[-1] + alpha[1])))
  trt = rbinom(N,1,ps)
  
  D = rbinom(N,1,ps)
  
  pY10 =  1/(1+exp(-drop(X %*% beta1[-1]+beta1[1])))
  pY00 = 1/(1+exp(- drop(X %*% beta0[-1]+beta0[1])))
  pY11 =  1/(1+exp(-drop(X %*% beta1[-1]-rev(beta1[1]))))
  pY01 = 1/(1+exp(- drop(X %*% beta0[-1]-rev(beta0[1]))))
  pY = trt*(1-D)*pY10 + (1-trt)*(1-D)*pY00 + trt*D*pY11 + (1-trt)*D*pY01
  Y = rbinom(N,1,pY)
  
  ate = mean((1-ps)*pY10 + ps*pY11 - (1-ps)*pY00 - ps*pY01)
  
  S = cbind( trt*st.mu + rnorm(N),
             Y*sy.mu + rnorm(N))
  
  return(list(Y = Y, trt = trt, 
              X = X, S = S,
              ate = ate))  
}

sim.demyst = function(N,p,alpha,beta1,beta0,
                        st.mu=1,sy.mu=1)
{
  X = matrix(rnorm(N*p),N,p)
  ps = 1/(1+exp(-drop(X %*% alpha[-1] + alpha[1])))
  trt = rbinom(N,1,ps)
  
  Z = cbind(exp(X[,1]/2),(X[,1]*X[,2]/25+0.6)^3)
  
  pY1 =  1/(1+exp(-drop(Z %*% beta1[-1]+beta1[1])))
  pY0 = 1/(1+exp(- drop(Z %*% beta0[-1]+beta0[1])))
  Y1 = rbinom(N,1,pY1)
  Y0 = rbinom(N,1,pY0)
  Y = trt*Y1 + (1-trt)*Y0
  
  S = cbind( trt*st.mu + rnorm(N),
             Y*sy.mu + rnorm(N))
  
  return(list(Y = Y, trt = trt, 
              X = X, S = S,
              ate = mean(Y1-Y0)))
}