require(splines)
source("core/utility.R")

cf.bs = function(Y, X, nfold, foldid, 
                 degree = 1, 
                 p.range = degree:round(sqrt(length(Y))^(1/ncol(X))),
                 fit.fun = function(formula,...){coef(lm(formula,...))}, 
                 subset = rep(T,length(Y)), 
                 dev.fun = ifelse(all(Y %in% 0:1), dev.binary.np, dev.ls),
                 link = identity, 
                 max.p.range = 10, 
                 ...)
{
  # Testing parameters
  # dat = sim.gen.np(1000)
  # Y = dat$Y[1:n]
  # X = dat$X
  # subset = dat$A[1:n]==1
  # p.range = 3:round(sqrt(length(Y)))
  # fit.fun = function(formula,...){coef(lm(formula))}
  # dev.fun = dev.ls
  # degree = 1
  
  if(length(p.range)>max.p.range)
  {
    p.range = round(seq(p.range[1],rev(p.range)[1], length.out = max.p.range) )
  }
  
  # body
  n = length(Y)
  N = nrow(X)
  sub.pos = which(subset)
  
  dev.best = dev.fun(Y[subset],rep(mean(Y[subset]),sum(subset)))
  p.best = 0
  full.pred = rep(mean(Y),N)
  cf.pred = rep(NA,N)
  for (ifold in 1:nfold)
  {
    cf.pred[foldid==ifold] = mean(Y[foldid[sub.pos]!=ifold])
  }
  
  bs.names = paste0("bs.X",formatC(1:ncol(X),width = 1+floor(log(ncol(X),10)), flag = 0))
  for(tensor in c(':','*'))
  {
    for (p in p.range)
    {
      for (j in 1:ncol(X))
      {
        assign(bs.names[j], bs(X[,j],degree=degree, df = p))
      }
      fit.formula = formula(paste0("~",paste0(bs.names, collapse = tensor)))
      bs.X = model.matrix(fit.formula)[,-1, drop = F]
      
      cf.coef = matrix(NA, ncol(bs.X)+1 ,nfold)
      cv.pred = rep(NA,n)
      for (ifold in 1:nfold)
      {
        trainpos = which(foldid[1:n]!= ifold & subset)
        cf.coef[,ifold] = fit.fun(formula(Y[trainpos]~bs.X[trainpos,]), ...)
        cvpos = which(foldid[1:n]== ifold)
        cf.coef[is.na(cf.coef[,ifold]),ifold] = 0
        cv.pred[cvpos] = link(cf.coef[1,ifold] + drop(bs.X[cvpos,, drop = F] %*% 
                                                   cf.coef[-1,ifold]))
      }
      tmp.dev = dev.fun(Y[subset],cv.pred[subset])
      
      if(tmp.dev < dev.best)
      {
        dev.best = tmp.dev
        p.best = p
        for (ifold in 1:nfold) 
        {
          cfpos = which(foldid== ifold)
          cf.pred[cfpos] = link(cf.coef[1,ifold] + drop(bs.X[cfpos,, drop = F] %*%
                                                     cf.coef[-1,ifold]))
        }
        full.coef = fit.fun(formula(Y[subset]~bs.X[sub.pos,]), ...)
        full.coef[is.na(full.coef)] = 0
        full.pred = link(full.coef[1]+drop(bs.X%*%full.coef[-1]))
      }
      # print(c(p,tmp.dev, p.best, dev.best))
      
    }}
  
  if(all(Y %in% 0:1))
  {
    meas = c(auc = auc(Y[subset], cf.pred[sub.pos],
                        levels = c(0,1), direction = "<"))
  }else{
    meas = c(Rsq = 1 - mean((Y[subset] - cf.pred[sub.pos])^2)/var(Y[subset]))
  }
  return(list(cf.pred = cf.pred, full.pred = full.pred, 
              df = p.best, dev = dev.best,
              meas = meas))
  
}

require(randomForest)

cf.rf = function(Y, X, nfold, foldid, subset = rep(T,length(Y)), 
                 dev.fun = ifelse(all(Y %in% 0:1), dev.binary.np, dev.ls),
                 ...)
{
  # Testing parameters
  Y = dat$Y[1:n]
  X = dat$X
  subset = dat$A[1:n]==1
  p.range = 3:round(sqrt(length(Y)))
  fit.fun = function(formula,...){coef(lm(formula))}
  dev.fun = dev.ls
  degree = 1
  
  # body
  n = length(Y)
  N = nrow(X)
  sub.pos = which(subset)
  if(all(class(X)!="data.frame"))
  {
    X = data.frame(X)
  }
  
  options(warn = -1) 
  # print(length(Y[subset]))
  # print(nrow(X[sub.pos,,drop=F]))
  full.fit = randomForest(x=X[sub.pos,,drop=F], y=Y[subset], xtest = X)
  options(warn = 1) 
  full.pred = full.fit$test$predicted
  cf.pred = rep(NA,N)
  for (ifold in 1:nfold)
  {
    trainpos = which(foldid[1:n]!=ifold & subset)
    testpos = which(foldid[1:n]==ifold & subset) 
    options(warn = -1) 
    # print(length(Y[trainpos]))
    # print(length(X[trainpos,]))
    cf.fit = randomForest(x=X[trainpos,,drop=F], y=Y[trainpos],
                          xtest = X[testpos,,drop=F],
                          ytest = Y[testpos],
                          keep.forest = T)
    options(warn = 1) 
    
    predpos = which(foldid==ifold) 
    cf.pred[predpos] = predict(cf.fit, X[predpos,,drop=F])
  }
  dev.rf = dev.fun(Y[subset], cf.pred[sub.pos])

  if(all(Y %in% 0:1))
  {
    meas = c(auc = auc(Y[subset], cf.pred[sub.pos],
                       levels = c(0,1), direction = "<"))
  }else{
    meas = c(Rsq = 1 - mean((Y[subset] - cf.pred[sub.pos])^2)/var(Y[subset]))
  }
  
  return(list(cf.pred = cf.pred, full.pred = full.pred, 
              dev = dev.rf,
              meas = meas))
}