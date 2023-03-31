SSL.aipw = function(Y,trt,X,S, XN,SN)
{
  n = length(Y)
  alpha = coef(glm(trt~X,family = binomial))
  beta1 = coef(glm(Y~X,family = binomial,subset = (trt==1))) 
  beta0 = coef(glm(Y~X,family = binomial,subset = (trt==0)))
  
  ps = expit(drop(alpha[1]+X %*% alpha[-1]))
  
  eta1 = coef(glm(Y~cbind(X,S,1/ps),family = binomial,subset = (trt==1))) 
  eta0 = coef(glm(Y~cbind(X,S,1/(1-ps)),family = binomial,subset = (trt==0))) 
  
  pred1 = expit(drop(beta1[1]+X%*% beta1[-1]))
  pred0 = expit(drop(beta0[1]+X%*% beta0[-1]))
  imp1 = expit(drop(eta1[1]+cbind(X,S,1/ps)%*% eta1[-1]))
  imp0 = expit(drop(eta0[1]+cbind(X,S,1/(1-ps))%*% eta0[-1]))
  Ut = (imp1-pred1)/ps + (imp0-pred0)/(1-ps)
  gr = coef(glm(trt~cbind(X,S,Ut),family = binomial))
  impt =  expit(drop(gr[1]+cbind(X,S,Ut)%*% gr[-1]))
  
  IX = cbind(1,X)
  S.alpha = (trt-ps)*IX
  S.beta1 = trt*(Y-pred1)*IX
  S.beta0 = (1-trt)*(Y-pred0)*IX
  S.gr = (trt-impt)*cbind(IX,S,Ut)
  S.eta1 = trt*(Y-imp1)*cbind(IX,S,1/ps)
  S.eta0 = (1-trt)*(Y-imp0)*cbind(IX,S,1/(1-ps))
  
  meat = var(cbind(S.alpha,S.beta1,S.beta0,S.gr,S.eta1,S.eta0))
  
  expps = exp(drop(alpha[1]+X %*% alpha[-1]))
  
  H.alpha = t(IX) %*% ( ps*(1-ps) * IX)/n
  H.beta1 = t(IX) %*% ( pred1*(1-pred1) * trt * IX)/n
  H.beta0 = t(IX) %*% ( pred0*(1-pred0) * (1-trt) * IX)/n
  H.gr.gr = t(cbind(IX,S,Ut)) %*% (impt*(1-impt)*cbind(IX,S,Ut))/n
  grUt = rev(gr)[1]
  H.gr.alpha = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-(imp1-pred1)/expps+(imp0-pred0)*expps)*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-(imp1-pred1)/expps+(imp0-pred0)*expps)*IX,
          2, mean)
  )
  H.gr.beta1 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-pred1*(1-pred1)/ps)*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-pred1*(1-pred1)/ps)*IX,
          2, mean)
  )
  H.gr.beta0 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-pred0*(1-pred0)/(1-ps))*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-pred0*(1-pred0)/(1-ps))*IX,
          2, mean)
  )
  H.gr.eta1 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (imp1*(1-imp1)/ps)*cbind(IX,S,1/ps))/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (imp1*(1-imp1)/ps)*cbind(IX,S,1/ps),
          2, mean)
  )
  H.gr.eta0 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (imp0*(1-imp0)/(1-ps))*cbind(IX,S,1/(1-ps)))/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (imp0*(1-imp0)/(1-ps))*cbind(IX,S,1/(1-ps)),
          2, mean)
  )
  H.eta1 = t(cbind(IX,S,1/ps)) %*% (trt*imp1*(1-imp1)*cbind(IX,S,1/ps))/n
  eta1u = rev(eta1)[1]
  H.eta1.alpha = rbind(
    t(cbind(IX,S)) %*% (trt*imp1*(1-imp1)*eta1u/(-expps)*IX)/n,
    apply(trt*(imp1*(1-imp1)/ps*eta1u - (Y-imp1))/(-expps)*IX,
          2,mean)
  )
  H.eta0 = t(cbind(IX,S,1/(1-ps))) %*% ((1-trt)*imp0*(1-imp0)*cbind(IX,S,1/(1-ps)))/n
  eta0u = rev(eta0)[1]
  H.eta0.alpha = rbind(
    t(cbind(IX,S)) %*% ((1-trt)*imp0*(1-imp0)*eta0u*expps*IX)/n,
    apply((1-trt)*(imp0*(1-imp0)/(1-ps)*eta0u - (Y-imp0))*expps*IX,
          2,mean)
  )
  
  p = ncol(IX)
  q = p + ncol(S) + 1
  negHess = matrix(0, p*3+q*3 , p*3+q*3)
  
  negHess[1:p,1:p] = H.alpha
  negHess[p+1:p,p+1:p] = H.beta1
  negHess[2*p+1:p,2*p+1:p] = H.beta0
  negHess[3*p+1:q,1:p] = H.gr.alpha
  negHess[3*p+1:q,p+1:p] = H.gr.beta1
  negHess[3*p+1:q,2*p+1:p] = H.gr.beta0
  negHess[3*p+1:q,3*p+1:q] = H.gr.gr
  negHess[3*p+1:q,3*p+q+1:q] = H.gr.eta1
  negHess[3*p+1:q,3*p+2*q+1:q] = H.gr.eta0
  negHess[3*p+q+1:q,1:p] = H.eta1.alpha
  negHess[3*p+q+1:q,3*p+q+1:q] = H.eta1
  negHess[3*p+2*q+1:q,1:p] = H.eta0.alpha
  negHess[3*p+2*q+1:q,3*p+2*q+1:q] = H.eta0
  
  ps = expit(drop(alpha[1]+XN %*% alpha[-1]))
  pred1 = expit(drop(beta1[1]+XN%*% beta1[-1]))
  pred0 = expit(drop(beta0[1]+XN%*% beta0[-1]))
  imp1 = expit(drop(eta1[1]+cbind(XN,SN,1/ps)%*% eta1[-1]))
  imp0 = expit(drop(eta0[1]+cbind(XN,SN,1/(1-ps))%*% eta0[-1]))
  Ut = (imp1-pred1)/ps + (imp0-pred0)/(1-ps) 
  impt = expit(drop(gr[1]+cbind(XN,SN,Ut)%*% gr[-1]))
  
  ate = mean(pred1 - pred0 + impt*(imp1-pred1)/ps - (1-impt)*(imp0-pred0)/(1-ps))
  
  IX = cbind(1,XN)
  ps1 = ps
  ps0 = 1-ps
  dps0 = exp(drop(alpha[1]+XN %*% alpha[-1]))
  dps1 = -1/dps0
  dpred1 = pred1*(1-pred1)
  dpred0 = pred0*(1-pred0)
  dimp1 = imp1*(1-imp1)
  dimp0 = imp0*(1-imp0)
  dimpt = impt*(1-impt)
  
  grad.alpha = apply((impt*(imp1-pred1)*dps1-(1-impt)*(imp0-pred0)*dps0
                      + impt*dimp1/ps1*eta1u*dps1 
                      - (1-impt)*dimp0/ps0*eta0u*dps0
                      + Ut*dimpt*grUt*((imp1-pred1)*dps1+(imp0-pred0)*dps0)
                    )*IX,2,mean)
  grad.beta1 = apply(((1-impt/ps1)*dpred1
                      + Ut*dimpt*grUt*(-dpred1/ps1))*IX,2,mean)
  grad.beta0 = apply((-(1-(1-impt)/ps0)*dpred0
                      + Ut*dimpt*grUt*(-dpred0/ps0))*IX,2,mean)
  grad.gr = apply(Ut*dimpt*cbind(IX,SN,Ut),2,mean)
  grad.eta1 = apply((impt/ps1*dimp1
                     + Ut*dimpt*grUt*(dimp1/ps1))*cbind(IX,SN,1/ps1),2,mean)
  grad.eta0 = apply((-(1-impt)/ps0*dimp0
                     + Ut*dimpt*grUt*(dimp0/ps0))*cbind(IX,SN,1/ps0),2,mean)
  
  grad = c(grad.alpha,grad.beta1,grad.beta0,grad.gr,grad.eta1,grad.eta0)
  u = solve(t(negHess),grad)
  
  se = drop(sqrt(t(u) %*% meat %*% u / n))
  
  return(list(est = ate, se =se))
}