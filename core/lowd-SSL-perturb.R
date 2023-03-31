SSL.aipw.perturb = function(Y,trt,X,S, XN,SN, nperturb)
{
  n = length(Y)
  N = nrow(XN)
  ate = rep(0, nperturb)
  for(iperturb in 1:nperturb)
  {
    G  = rbeta(n,0.5,1.5)*4
    GN  = rbeta(N,0.5,1.5)*4
    alpha = coef(glm(trt~X,family = binomial,weights = G))
    beta1 = coef(glm(Y~X,family = binomial,subset = (trt==1),weights = G)) 
    beta0 = coef(glm(Y~X,family = binomial,subset = (trt==0),weights = G))
    
    ps = expit(drop(alpha[1]+X %*% alpha[-1]))
    
    eta1 = coef(glm(Y~cbind(X,S,1/ps),family = binomial,subset = (trt==1))) 
    eta0 = coef(glm(Y~cbind(X,S,1/(1-ps)),family = binomial,subset = (trt==0))) 
    
    pred1 = expit(drop(beta1[1]+X%*% beta1[-1]))
    pred0 = expit(drop(beta0[1]+X%*% beta0[-1]))
    imp1 = expit(drop(eta1[1]+cbind(X,S,1/ps)%*% eta1[-1]))
    imp0 = expit(drop(eta0[1]+cbind(X,S,1/(1-ps))%*% eta0[-1]))
    Ut = (imp1-pred1)/ps + (imp0-pred0)/(1-ps)
    gr = coef(glm(trt~cbind(X,S,Ut),family = binomial, weights = G))
    
    ps = expit(drop(alpha[1]+XN %*% alpha[-1]))
    pred1 = expit(drop(beta1[1]+XN%*% beta1[-1]))
    pred0 = expit(drop(beta0[1]+XN%*% beta0[-1]))
    imp1 = expit(drop(eta1[1]+cbind(XN,SN,1/ps)%*% eta1[-1]))
    imp0 = expit(drop(eta0[1]+cbind(XN,SN,1/(1-ps))%*% eta0[-1]))
    Ut = (imp1-pred1)/ps + (imp0-pred0)/(1-ps) 
    impt = expit(drop(gr[1]+cbind(XN,SN,Ut)%*% gr[-1]))
    
    ate[iperturb] = mean((pred1 - pred0 + impt*(imp1-pred1)/ps - (1-impt)*(imp0-pred0)/(1-ps))
              * GN)
  }
  return(list(se=sd(ate), CI = quantile(ate,probs = c(.025,.975))))
}