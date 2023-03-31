
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
  non.NA = !apply(is.na(full.fit$bet.all), 2, any)
  full.fit$rho = full.fit$rho[non.NA]
  full.fit$bet.all = full.fit$bet.all[,non.NA]
  
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
  
  full.coef = full.fit$bet.all[,best.ilam]
  full.pred = expit(drop(X %*% full.coef[-1])+full.coef[1])
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
              full.coef = full.coef, 
              lambda = full.fit$rho[best.ilam], dev = min(apply(dev.cv,2,mean)),
              meas = meas))
  
}