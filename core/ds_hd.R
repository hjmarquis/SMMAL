
source("core/utility.R")

require(glmnet)
ds.lasso = function(Y, X, nfold, foldid, dsgrp,
                 subset = rep(T,length(Y)), 
                 family = ifelse(all(Y %in% 0:1), "binomial", "gaussian"),
                 dev.fun = ifelse(all(Y %in% 0:1), dev.logistic, dev.ls),
                 weight = rep(1,nrow(X)), 
                 ds.weight = NULL,
                 ...)
{
  
  n = length(Y)
  N = nrow(X)
  if(hasArg("penalty.factor"))
  {
    if(all(is.infinite(list(...)$penalty.factor)))
    {
      return(cf.null(Y, ncol(X), nfold, foldid, subset, family,weight[1:n]))
    }
  }
  
  
  # body
  link = ifelse(family == "binomial", expit, identity)

  sub.pos = which(subset & dsgrp[1:n])
  full.fit = glmnet(X[sub.pos,], Y[sub.pos], family = family,
                    weights = weight[sub.pos],...)
  cv.fit = vector("list", nfold)
  
  ifold = 1
  dev.cv = matrix(NA, nfold, length(full.fit$lambda))
  for(ifold in 1:nfold)
  {
    train.pos = which(subset & (foldid[1:n]!=ifold) & dsgrp[1:n])
    if(!is.null(ds.weight))
    {
      train.wgt = ds.weight[train.pos,ifold]
    }else{
      train.wgt = weight[train.pos]
    }
    cv.fit[[ifold]] = glmnet(X[train.pos,], Y[train.pos], family = family,
                             lambda = full.fit$lambda,
                             weights = train.wgt,...)
    test.pos = which(subset & (foldid[1:n]==ifold) & dsgrp[1:n])
    dev.cv[ifold,] = dev.fun(Y[test.pos], X[test.pos,] %*% cv.fit[[ifold]]$beta
                             + rep(cv.fit[[ifold]]$a0, each = length(test.pos)),
                             weight[test.pos])
  }
  best.lam = full.fit$lambda[which.min(apply(dev.cv,2,mean))]

  tmp.coef = drop(coef(full.fit, s = best.lam))
  full.pred = link(drop(X %*% tmp.coef[-1])+tmp.coef[1])
  cf.pred = rep(0, N)
  ds.pred = matrix(NA,N,nfold)
  
  for (ifold in 1:nfold)
  {
    test.pos = which(foldid==ifold)
    tmp.coef = drop(coef(cv.fit[[ifold]], s = best.lam))
    cf.pred[test.pos] = link(drop(X[test.pos,] %*% tmp.coef[-1])+
                               tmp.coef[1])
    ds.pos = which(foldid!=ifold & !dsgrp)
    ds.pred[ds.pos, ifold] = link(drop(X[ds.pos,] %*% tmp.coef[-1])+
                                    tmp.coef[1])
  }
  
  if(all(Y %in% 0:1))
  {
    meas = c(auc = auc(Y[sub.pos], cf.pred[sub.pos],
                        levels = c(0,1), direction = "<"))
  }else{
    meas = c(Rsq = 1 - mean((Y[sub.pos] - cf.pred[sub.pos])^2)/var(Y[sub.pos]))
  }
  
  out = list(cf.pred = cf.pred, full.pred = full.pred, 
             full.coef = tmp.coef, 
             lambda = best.lam, dev = min(apply(dev.cv,2,mean)),
             meas = meas, ds.pred = ds.pred)
  
  return(out)
  
}
