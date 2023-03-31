
  n = length(Y)
# Initial estimators
  alpha = coef(cv.glmnet(X,trt, family="binomial"), s = 'lambda.min')
  beta1 = coef(cv.glmnet(X[trt==1,],Y[trt==1], family="binomial"), s = 'lambda.min')
  beta0 = coef(cv.glmnet(X[trt==0,],Y[trt==0], family="binomial"), s = 'lambda.min')
  
  ps = expit(drop(alpha[1]+X %*% alpha[-1]))

# Imputation for estimation  
  eta1 = coef(cv.glmnet(cbind(X,S,1/ps)[trt==1,],Y[trt==1], family="binomial"), s = 'lambda.min')
  eta0 = coef(cv.glmnet(cbind(X,S,1/(1-ps))[trt==0,],Y[trt==0], family="binomial"), s = 'lambda.min')
  
  pred1 = expit(drop(beta1[1]+X%*% beta1[-1]))
  pred0 = expit(drop(beta0[1]+X%*% beta0[-1]))
  imp1 = expit(drop(eta1[1]+cbind(X,S,1/ps)%*% eta1[-1]))
  imp0 = expit(drop(eta0[1]+cbind(X,S,1/(1-ps))%*% eta0[-1]))
  Ut = (imp1-pred1)/ps + (imp0-pred0)/(1-ps)
  gr = coef(glm(trt~cbind(X,S,Ut),family = binomial))
  impt =  expit(drop(gr[1]+cbind(X,S,Ut)%*% gr[-1]))

# Estimation  
  ate = mean(pred1 - pred0 + impt*(imp1-pred1)/ps - (1-impt)*(imp0-pred0)/(1-ps))
  
# Imputation for Hessian
  
# Hessian
  
  expps = exp(drop(alpha[1]+X %*% alpha[-1]))
  
  H.alpha = t(IX) %*% ( ps*(1-ps) * IX)/n
  # with imputed trt ~ pred1*(1-pred1)*IX^2
  H.beta1 = t(IX) %*% ( pred1*(1-pred1) * trt * IX)/n
  # with imputed trt ~ pred0*(1-pred0)*IX^2
  H.beta0 = t(IX) %*% ( pred0*(1-pred0) * (1-trt) * IX)/n
  H.gr.gr = t(cbind(IX,S,Ut)) %*% (impt*(1-impt)*cbind(IX,S,Ut))/n
  grUt = rev(gr)[1]
  # with imputed trt~ (-(imp1-pred1)/expps+(imp0-pred0)*expps)*IX
  H.gr.alpha = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-(imp1-pred1)/expps+(imp0-pred0)*expps)*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-(imp1-pred1)/expps+(imp0-pred0)*expps)*IX,
          2, mean)
  )
  #  with imputed trt~ (-pred1*(1-pred1)/ps)*IX
  H.gr.beta1 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-pred1*(1-pred1)/ps)*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-pred1*(1-pred1)/ps)*IX,
          2, mean)
  )
  #  with imputed trt~ (-pred0*(1-pred0)/(1-ps))*IX
  H.gr.beta0 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (-pred0*(1-pred0)/(1-ps))*IX)/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (-pred0*(1-pred0)/(1-ps))*IX,
          2, mean)
  )
  # with imputed trt~(imp1*(1-imp1)/ps)*cbind(IX,S,1/ps)
  H.gr.eta1 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (imp1*(1-imp1)/ps)*cbind(IX,S,1/ps))/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (imp1*(1-imp1)/ps)*cbind(IX,S,1/ps),
          2, mean)
  )
  # with imputed trt~(imp0*(1-imp0)/(1-ps))*cbind(IX,S,1/(1-ps))
  H.gr.eta0 = rbind(
    t(cbind(IX,S)) %*% (impt*(1-impt)*grUt*
                          (imp0*(1-imp0)/(1-ps))*cbind(IX,S,1/(1-ps)))/n, 
    apply((Ut*impt*(1-impt)*grUt-(trt-impt))*
            (imp0*(1-imp0)/(1-ps))*cbind(IX,S,1/(1-ps)),
          2, mean)
  )
  
  # with imputed trt~imp1*(1-imp1)*cbind(IX,S,1/ps)^2
  H.eta1 = t(cbind(IX,S,1/ps)) %*% (trt*imp1*(1-imp1)*cbind(IX,S,1/ps))/n
  eta1u = rev(eta1)[1]
  # with imputed trt~imp1*(1-imp1)*eta1u/(-expps)*IX*(IX,S)
  # with imputed Y ~ trt/(-expps)*IX
  # with imputed trt~(imp1*(1-imp1)/ps*eta1u - (Y-imp1))/(-expps)*IX
  H.eta1.alpha = rbind(
    t(cbind(IX,S)) %*% (trt*imp1*(1-imp1)*eta1u/(-expps)*IX)/n,
    apply(trt*(imp1*(1-imp1)/ps*eta1u - (Y-imp1))/(-expps)*IX,
          2,mean)
  )
  # with imputed trt~imp0*(1-imp0)*cbind(IX,S,1/(1-ps))^2
  H.eta0 = t(cbind(IX,S,1/(1-ps))) %*% ((1-trt)*imp0*(1-imp0)*cbind(IX,S,1/(1-ps)))/n
  eta0u = rev(eta0)[1]
  # with imputed trt~imp0*(1-imp0)*eta0u*expps*IX
  # with imputed Y ~ (1-trt)*expps*IX
  # with imputed trt ~ (imp0*(1-imp0)/(1-ps)*eta0u - (Y-imp0))*expps*IX
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
  
# Bias correction
  
# Variance estimation
