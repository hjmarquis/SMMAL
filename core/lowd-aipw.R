aipw = function(Y,trt,X)
{
  n = length(Y)
  alpha = coef(glm(trt~X,family = binomial))
  beta1 = coef(glm(Y~X,family = binomial,subset = (trt==1))) 
  beta0 = coef(glm(Y~X,family = binomial,subset = (trt==0)))
  
  pred1 = expit(beta1[1]+drop(X%*%beta1[-1]))
  pred0 = expit(beta0[1]+drop(X%*%beta0[-1]))
  ps1 = expit(alpha[1]+drop(X %*% alpha[-1]))
  ps0 = 1-ps1
  
  aipw = pred1 - pred0 + trt*(Y-pred1)/ps1 - (1-trt)*(Y-pred0)/ps0
  
  return(list(est = mean(aipw),
              se = sd(aipw)/sqrt(n)))
}

aipw.oracle = function(Y,trt,X, alpha, beta1, beta0)
{
  n = length(Y)
  
  pred1 = expit(beta1[1]+drop(X%*%beta1[-1]))
  pred0 = expit(beta0[1]+drop(X%*%beta0[-1]))
  ps1 = expit(alpha[1]+drop(X %*% alpha[-1]))
  ps0 = 1-ps1
  
  aipw = pred1 - pred0 + trt*(Y-pred1)/ps1 - (1-trt)*(Y-pred0)/ps0
  
  return(list(est = mean(aipw),
              se = sd(aipw)/sqrt(n)))
}
