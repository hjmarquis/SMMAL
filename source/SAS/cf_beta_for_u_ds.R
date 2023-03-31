# Y = data$Y[1:n]
# X = data$X
# S = data$S


cf.beta.u.ds = function(Y,X,S,gamma,nfold,nfoldid,Nfoldid,
                        penalty.factor = rep(1,ncol(X)),
                        u.ds = TRUE, nopen = FALSE)
{
  n = length(Y)
  p = ncol(X)
  N = nrow(X) - n
  W = cbind(X,S)
  min.pen = ifelse((2*p)<(N+n), 1e-7, 1e-3)
  
  Ydag.full = c(expit(gamma$full[1]+
                 drop(W[-1:-n,]%*%gamma$full[-1])))
  # start.time = Sys.time()
  full.nopen = glm(c(Y,Ydag.full)~X, family = binomial)
  # print(Sys.time()-start.time)
  
  dev.nopen.fold = rep(0,nfold)
  cf.nopen.fit = vector("list",nfold)
  Ydag.cf = c(Y, rep(0, N))
  for (ifold in 1:nfold)
  {
    train.npos = which(nfoldid != ifold)
    train.Npos = n+which(Nfoldid != ifold)
    Ydag = c(Y[train.npos],expit(gamma$cf[1,ifold]+
                   drop(W[train.Npos,]%*%gamma$cf[-1,ifold])))
    # start.time = Sys.time()
    cf.nopen.fit[[ifold]] = glm(Ydag~X[c(train.npos, train.Npos),],family = binomial,
                             start = coef(full.nopen))
    # print(Sys.time()-start.time)
    
    test.npos = which(nfoldid == ifold)
    test.Npos = n+which(Nfoldid == ifold)
    pred = cbind(1,X[c(test.npos,test.Npos),]) %*% 
      coef(cf.nopen.fit[[ifold]])
    Ydag.cf[test.Npos] = expit(gamma$cf[1,ifold]+
                                 drop(W[test.Npos,]%*%gamma$cf[-1,ifold]))
    dev.nopen.fold[ifold] = logistic.dev(Ydag.cf[c(test.npos,test.Npos)],pred)
  }
  full = coef(full.nopen)
  cf = sapply(cf.nopen.fit,coef)
  best.lambda = 0
  
  if(!nopen)
  {
    Xdummy = X[c(1:n,rep(n+1:N,2)),]
    Ydummy = c(Y, rep(1:0, each = N))
    wdummy = c(rep(1,n), Ydag.full, 1-Ydag.full)
    lambda.list = exp(seq(log(glmnet(Xdummy,Ydummy,family = "binomial",
                         weights = wdummy, nlambda = 3)$lambda[1]),
                      log(min.pen), length.out = 100))
    full.fit = glmnet(Xdummy,Ydummy,family = "binomial",
                         weights = wdummy, lambda = lambda.list,
                      penalty.factor = penalty.factor)
    
    cf.fit = vector("list",nfold)
    dev.fold = matrix(0,nfold,length(full.fit$lambda))
    
    for (ifold in 1:nfold)
    {
      train.npos = which(nfoldid != ifold)
      train.Npos = n+which(Nfoldid != ifold)
      Ydag = expit(gamma$cf[1,ifold]+
                     drop(W[train.Npos,]%*%gamma$cf[-1,ifold]))
      Xdummy = X[c(train.npos,rep(train.Npos,2)),]
      Ydummy = c(Y[train.npos], rep(1:0, each = length(train.Npos)))
      wdummy = c(rep(1,length(train.npos)), Ydag, 1-Ydag)
      cf.fit[[ifold]] = glmnet(Xdummy,Ydummy,family = "binomial",
                               weights = wdummy, lambda = lambda.list,
                               penalty.factor = penalty.factor)
      
      test.npos = which(nfoldid == ifold)
      test.Npos = n+which(Nfoldid == ifold)
      pred = cbind(1,X[c(test.npos,test.Npos),]) %*% 
                     coef(cf.fit[[ifold]],s=full.fit$lambda)
      Ydag.cf[test.Npos] = expit(gamma$cf[1,ifold]+
                     drop(W[test.Npos,]%*%gamma$cf[-1,ifold]))
      dev.fold[ifold,] = logistic.dev(Ydag.cf[c(test.npos,test.Npos)],pred)
    }
    best.dev = min(apply(dev.fold,2,mean))
    
    if(best.dev < mean(dev.nopen.fold))
    {
      best.lambda = full.fit$lambda[which.min(apply(dev.fold,2,mean))]/sqrt(nfold/(nfold-1))
      full = as.numeric(coef(full.fit,s=best.lambda))
      cf = sapply(cf.fit,coef, s = best.lambda*sqrt(nfold/(nfold-1)))
      cf = sapply(cf,drop)
    }
  }
  if(!u.ds)
  {
    return(list(full=full,cf=cf,lambda = best.lambda))
  }
  
  res.cf = rep(NA, N+n)
  wgt.cf = rep(NA, N+n)
  wgt.cf.cf = matrix(NA,N+n,nfold)
  
  for (ifold in 1:nfold)
  {
    test.pos = c(which(nfoldid == ifold),
                n+which(Nfoldid == ifold))
    
    pred = expit(cf[1,ifold] + drop(X[test.pos,] %*% cf[-1,ifold]))
    res.cf[test.pos] = pred - Ydag.cf[test.pos]
    wgt.cf[test.pos] = pred*(1-pred)
    
    for(jfold in setdiff(1:nfold, ifold))
    {
      train.npos = which(!(nfoldid %in% c(ifold, jfold)))
      train.Npos = n+which(!(Nfoldid %in% c(ifold, jfold)))
      tmpfit = glmnet(W[train.npos,], Y[train.npos], family = "binomial")
      tmpcoef = as.numeric(coef(tmpfit,s=gamma$lambda*sqrt(nfold/(nfold-2))))
      Ydag = expit(tmpcoef[1]+
                     drop(W[train.Npos,]%*%tmpcoef[-1]))
      if(best.lambda == 0)
      {
        tmpfit = glm(c(Y[train.npos], Ydag)~X[c(train.npos, train.Npos),],
                     family = binomial, start = coef(full.nopen))
        tmpcoef = coef(tmpfit)
      }else{
        Xdummy = X[c(train.npos,rep(train.Npos,2)),]
        Ydummy = c(Y[train.npos], rep(1:0, each = length(train.Npos)))
        wdummy = c(rep(1,length(train.npos)), Ydag, 1-Ydag)
        tmpfit = glmnet(Xdummy,Ydummy,family = "binomial",
                        weights = wdummy,
                        penalty.factor = penalty.factor)
        tmpcoef = as.numeric(coef(tmpfit,s=best.lambda*sqrt(nfold/(nfold-2))))
      }
      
      test.pos = c(which(nfoldid == jfold),
                    n+which(Nfoldid == jfold))
      pred = expit(tmpcoef[1] + drop(X[test.pos,] %*% tmpcoef[-1]))
      wgt.cf.cf[test.pos,ifold] = pred*(1-pred)
    }
  }
  
  return(list(full=full,cf=cf, wgt.cf = wgt.cf, wgt.cf.cf = wgt.cf.cf, 
              res.cf = res.cf, lambda = best.lambda))
}