# load("other/test/test_dat.rda")
# eta = eta$full
# 
# rm(list = setdiff(objects(), 
#                   c("X","S","trt","eta","n")))
# source("source/utility.R")
# maxit = 1000
# tol = 1e-7
# link = expit
# dlink = dexpit
# inv.link = logit
# nlambda = 100
# max.df = 100
# foldid = nfoldid
# max.step = 1
# init.xi1 = 0.5
# eta = eta1
# 
# eta = eta0
# trt = 1-trt
# penalty.factor = rep(c(1,0,1),c(300,1,11))

cv.beta.xi = function(X,S,trt,eta,n,
                      link, dlink, inv.link,
                      nfold, foldid, 
                      nlambda = 100, max.df = sqrt(n), 
                      maxit=1000, tol = 1e-7,
                      max.step = 1, 
                      init.xi1 = 0.5)
{
  full.fit = beta.xi.full(X,S,trt,eta,n,
               link, dlink, inv.link,
               nlambda = nlambda, max.df = max.df, 
               maxit=maxit, tol = tol, max.step = max.step,
               init.xi1 = init.xi1)
 
  N = nrow(X)
  cf.fit = vector("list",nfold)
  dev.fold = matrix(NA,nfold,length(full.fit$lambda))
  for(ifold in 1:nfold)
  {
    fold.pos = which(foldid!=ifold)
    fold.nlab = length(fold.pos)
    SSL.pos = c(fold.pos, (n+1):N)
    
    cf.fit[[ifold]] = beta.xi.full(X[SSL.pos,],S[SSL.pos,],trt[fold.pos],
                               eta,fold.nlab,
                               link, dlink, inv.link,
                               max.df = max.df, 
                               lambda.list = full.fit$lambda, 
                               maxit=maxit, tol = tol,
                               max.step = max.step,
                               beta = full.fit$beta[,1],
                               xi = full.fit$xi[,1])
    test.pos = which(foldid==ifold)
    cv.pred = link(cbind(1,X[test.pos,],S[test.pos,]) %*% cf.fit[[ifold]]$xi)
      
    dev.fold[ifold,1:length(cf.fit[[ifold]]$lambda)] = logistic.dev(trt[test.pos],cv.pred)
  }
  best.ilam = which.min(apply(dev.fold,2,mean))
  full = list(xi = full.fit$xi[,best.ilam],
              beta = full.fit$beta[,best.ilam],
              lambda = full.fit$lambda[best.ilam])
  cf = list(xi = sapply(cf.fit, function(x) x$xi[,best.ilam]),
            beta = sapply(cf.fit, function(x) x$beta[,best.ilam]),
            lambda = full.fit$lambda[best.ilam])
  return(list(full=full,cf=cf))
}