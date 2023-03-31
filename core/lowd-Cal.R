cal.aipw = function(Y,trt,X)
{
  n = length(Y)
  init.alpha = coef(glm(trt~X,family = binomial))
  init.beta1 = coef(glm(Y~X,family = binomial,subset = (trt==1))) 
  init.beta0 = coef(glm(Y~X,family = binomial,subset = (trt==0)))
  
  
  w.wd1 = trt*exp(-init.alpha[1]-drop(X %*% init.alpha[-1]))
  w.wd0 = (1-trt)*exp(init.alpha[1]+drop(X %*% init.alpha[-1]))

  wd.beta1 = coef(glm(Y~X,family = binomial, weights = w.wd1))
  wd.beta0 = coef(glm(Y~X,family = binomial, weights = w.wd0))

  init.pred1 = expit(init.beta1[1]+drop(X %*% init.beta1[-1]))
  w.cal1 = init.pred1*(1-init.pred1)
  init.pred0 = expit(init.beta0[1]+drop(X %*% init.beta0[-1]))
  w.cal0 = init.pred0*(1-init.pred0)
  
  cal1.alpha = cal.fit(trt,X,w.cal1,init.alpha)
  cal0.alpha = -cal.fit(1-trt,X,w.cal0,-init.alpha)
  
  pred1 = expit(wd.beta1[1]+drop(X%*%wd.beta1[-1]))
  pred0 = expit(wd.beta0[1]+drop(X%*%wd.beta0[-1]))
  ps1 = expit(cal1.alpha[1]+drop(X %*% cal1.alpha[-1]))
  ps0 = 1-expit(cal0.alpha[1]+drop(X %*% cal0.alpha[-1]))
  
  aipw = pred1 - pred0 + trt*(Y-pred1)/ps1 - (1-trt)*(Y-pred0)/ps0
  
  return(list(est = mean(aipw),
              se = sd(aipw)/sqrt(n)))
}