
source("core/utility.R")

cf.null = function(Y, p, nfold, foldid, 
                   subset = rep(T,length(Y)),
                   family = ifelse(all(Y %in% 0:1), "binomial", "gaussian"),
                   weight = rep(1,length(Y)))
{
  family = get(family)
  
  full.fit = glm(Y~1, family = family,
                 subset = subset, 
                    weights = weight)
  cf.pred = rep(0, length(foldid))
  
  for(ifold in 1:nfold)
  {
    cf.pred[foldid==ifold] = mean(Y[foldid[1:n]!=ifold])
  }
  
  return(list(cf.pred = cf.pred, 
              full.pred = rep(mean(Y),n), 
              full.coef = c(coef(full.fit),rep(0,p))))
}

require(glmnet)
cf.lasso = function(Y, X, nfold, foldid, 
                 subset = rep(T,length(Y)), 
                 family = ifelse(all(Y %in% 0:1), "binomial", "gaussian"),
                 dev.fun = ifelse(all(Y %in% 0:1), dev.logistic, dev.ls),
                 weight = rep(1,nrow(X)), 
                 dcf.weight = NULL, 
                 dcf = FALSE,
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

  sub.pos = which(subset)
  
  full.fit = glmnet(X[sub.pos,], Y[sub.pos], family = family,
                    weights = weight[sub.pos],...)
  cv.fit = vector("list", nfold)
  
  ifold = 1
  dev.cv = matrix(NA, nfold, length(full.fit$lambda))
  for(ifold in 1:nfold)
  {
    train.pos = which(subset & (foldid[1:n]!=ifold))
    if(!is.null(dcf.weight))
    {
      train.wgt = dcf.weight[train.pos,ifold]
    }else{
      train.wgt = weight[train.pos]
    }
    cv.fit[[ifold]] = glmnet(X[train.pos,], Y[train.pos], family = family,
                             lambda = full.fit$lambda,
                             weights = train.wgt,...)
    test.pos = which(subset & (foldid[1:n]==ifold))
    tmp.dev = dev.fun(Y[test.pos], X[test.pos,] %*% cv.fit[[ifold]]$beta
                             + rep(cv.fit[[ifold]]$a0, each = length(test.pos)),
                             weight[test.pos])
    dev.cv[ifold,1:length(tmp.dev)] = tmp.dev
  }
  best.lam = full.fit$lambda[which.min(apply(dev.cv,2,mean))]

  tmp.coef = drop(coef(full.fit, s = best.lam))
  full.pred = link(drop(X %*% tmp.coef[-1])+tmp.coef[1])
  cf.pred = rep(0, N)
  
  for (ifold in 1:nfold)
  {
    test.pos = which(foldid==ifold)
    tmp.coef = drop(coef(cv.fit[[ifold]], s = best.lam))
    cf.pred[test.pos] = link(drop(X[test.pos,] %*% tmp.coef[-1])+
                               tmp.coef[1])
  }
  
  if(all(Y %in% 0:1))
  {
    meas = c(auc = auc(Y[subset], cf.pred[sub.pos],
                        levels = c(0,1), direction = "<"))
  }else{
    meas = c(Rsq = 1 - mean((Y[subset] - cf.pred[sub.pos])^2)/var(Y[subset]))
  }
  
  out = list(cf.pred = cf.pred, full.pred = full.pred, 
             full.coef = tmp.coef, 
             lambda = best.lam, dev = min(apply(dev.cv,2,mean)),
             meas = meas)
  
  if(dcf)
  {
    out$dcf.pred = matrix(NA, N, nfold)
    for(ifold in 1:(nfold-1))
    {
      for (jfold in (ifold+1):nfold)
      {
        train.pos = which(subset & !(foldid[1:n] %in% c(ifold,jfold)))
        dcf.fit = glmnet(X[train.pos,], Y[train.pos], family = family,
                                 lambda = full.fit$lambda,
                                 weights = weight[train.pos],...)
        tmp.coef = drop(coef(dcf.fit, s = best.lam))
        
        test.pos = which(foldid==ifold)
        out$dcf.pred[test.pos, jfold] = link(drop(X[test.pos,] %*% tmp.coef[-1])+
                                   tmp.coef[1])
        test.pos = which(foldid==jfold)
        out$dcf.pred[test.pos, ifold] = link(drop(X[test.pos,] %*% tmp.coef[-1])+
                                           tmp.coef[1])
      }
    }
  }
  
  return(out)
  
}


require(RCAL)
# 
cf.rcal = function(Y, X, nfold, foldid,
                   subset = rep(T,length(Y)), 
                   weight = rep(1,nrow(X)),
                   dcf.weight = NULL,
                   ...)
{
  # Testing parameters
  # coef = c(0,1,rep(0,49))
  # dat = sim.gen.hd(1000, 50,
  #                  par.lp.A = list(coef=coef),
  #                  par.lp.Y1 = list(coef=coef),
  #                  par.lp.Y0 = list(coef=coef))
  # Y = dat$Y[1:300]
  # X = dat$X
  # nfold = 5
  # foldid = c(rep_len(1:5,300)[sample(300)],
  #            rep_len(1:5,700)[sample(700)])
  # subset = dat$A[1:300]==1
  # weight = rep(1,nrow(X))
  
  
  
  # body
  n = length(Y)
  N = nrow(X)
  sub.pos = which(subset)
  full.fit = glm.regu.path(nrho = 100, tune.fac = sqrt(0.95), 
                           x=X[sub.pos,], y=Y[sub.pos], 
                           iw = weight[sub.pos],...)
  cv.fit = vector("list", nfold)
  non.NA = which(!apply(is.na(full.fit$bet.all), 2, any))
  full.fit$rho = c(full.fit$rho[non.NA],
                   max(full.fit$rho,na.rm=TRUE)*2)
  full.fit$bet.all = full.fit$bet.all[
    , c(non.NA,max(non.NA))]
  
  ifold = 1
  dev.cv = matrix(NA, nfold, length(full.fit$rho))
  for(ifold in 1:nfold)
  {
    train.pos = which(subset & (foldid[1:n]!=ifold))
    if(!is.null(dcf.weight))
    {
      train.wgt = dcf.weight[train.pos,ifold]
    }else{
      train.wgt = weight[train.pos]
    }
    cv.fit[[ifold]] = glm.regu.path(rho.seq = full.fit$rho, 
                                    x=X[train.pos,], y=Y[train.pos], 
                                    iw = train.wgt, ...)
    test.pos = which(subset & (foldid[1:n]==ifold))
    dev.cv[ifold,] = dev.rcal(Y[test.pos], X[test.pos,] %*% cv.fit[[ifold]]$bet.all[-1,]
                             + rep(cv.fit[[ifold]]$bet.all[1,], each = length(test.pos)),
                             weight[test.pos])
  }
  mdev = apply(dev.cv,2,mean)
  if(all(is.na(mdev)))
  {
    mdev = apply(dev.cv,2,mean,na.rm = TRUE)
  }
  best.ilam = which.min(mdev)
  
  tmp.coef = full.fit$bet.all[,best.ilam]
  full.pred = expit(drop(X %*% tmp.coef[-1])+tmp.coef[1])
  cf.pred = rep(0, N)
  
  for (ifold in 1:nfold)
  {
    test.pos = which(foldid==ifold)
    tmp.coef = cv.fit[[ifold]]$bet.all[,best.ilam]
    cf.pred[test.pos] = expit(drop(X[test.pos,] %*% tmp.coef[-1])+
                               tmp.coef[1])
  }
  
  if(all(Y %in% 0:1))
  {
    meas = c(auc = auc(Y[subset], cf.pred[sub.pos],
                       levels = c(0,1), direction = "<"))
  }else{
    meas = c(Rsq = 1 - mean((Y[subset] - cf.pred[sub.pos])^2)/var(Y[subset]))
  }
  return(list(cf.pred = cf.pred, full.pred = full.pred, 
              full.coef = tmp.coef, 
              lambda = full.fit$rho[best.ilam], dev = min(apply(dev.cv,2,mean)),
              meas = meas))
  
}