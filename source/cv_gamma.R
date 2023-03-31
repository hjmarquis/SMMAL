cv.gamma = function(Y,X,S,nfold,foldid, weights = rep(1,length(Y)))
{
  n = length(Y)
  W = cbind(X,S)
  
  gfit = glmnet(W,Y,family = "binomial",weights = weights)
  gfit.cf = vector("list",nfold)
  dev.fold = matrix(0,nfold,length(gfit$lambda))
  for(ifold in 1:nfold)
  {
    gfit.cf[[ifold]] = glmnet(W[foldid!=ifold,],Y[foldid!=ifold],family = "binomial",
                              weights = weights[foldid!=ifold])
    cv.pred = predict(gfit.cf[[ifold]], newx = W[foldid==ifold,],
                      type = "response", s = gfit$lambda)
    dev.fold[ifold,] = logistic.dev(Y[foldid==ifold],cv.pred,
                                    weights = weights[foldid==ifold])
  }
  best.lambda = gfit$lambda[which.min(apply(dev.fold,2,mean))]
  full = coef(gfit,s=best.lambda)
  cf = sapply(gfit.cf,coef, s = best.lambda*sqrt(nfold/(nfold-1)))
  cf = sapply(cf,drop)
  
  return(list(full=full,cf=cf))
}