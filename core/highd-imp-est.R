highd.imp.est = function(Y,trt,X,S, init)
{
  ps = expit(drop(init$alpha[1]+X %*% init$alpha[-1]))
  
  eta1 = coef(cv.glmnet(cbind(X,S,1/ps)[trt==1,],Y[trt==1], family="binomial"), s = 'lambda.min')
  eta0 = coef(cv.glmnet(cbind(X,S,1/(1-ps))[trt==0,],Y[trt==0], family="binomial"), s = 'lambda.min')

  gr = coef(cv.glmnet(cbind(X,S,Ut), trt, family="binomial"), s = 'lambda.min')
  
  return(list(eta1=eta1,eta0=eta0,gr=gr))
}