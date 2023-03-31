highd.init =function(Y,trt,X,S)
{
  alpha = coef(cv.glmnet(X,trt, family="binomial"), s = 'lambda.min')
  beta1 = coef(cv.glmnet(X[trt==1,],Y[trt==1], family="binomial"), s = 'lambda.min')
  beta0 = coef(cv.glmnet(X[trt==0,],Y[trt==0], family="binomial"), s = 'lambda.min')
  
  return(list(alpha = alpha, 
         beta1 = beta1,
         beta0 = beta0))
}