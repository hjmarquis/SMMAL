sim.cov.emr2 = function(N,p,q,X.factor = 5, 
                         kterm = 100)
{
  U = rnorm(N)/sqrt(p)
  Xraw = matrix(rnorm(N*(p)),N,p)*sqrt(1-1/(p))+U 
  # X1raw = Xraw %*% rep(0.5/sqrt(p-1),p-1) + rnorm(N)*sqrt(.5)
  # Xraw = cbind(X1raw,Xraw)
  Sraw = matrix(rnorm(N*q),N,q)
  
  prob.X = diff(pnorm(log(exp(0:(kterm+1))-1)/X.factor))
  Xm = sum((0:kterm)*prob.X)
  X.sd = sqrt(sum((0:kterm)^2*prob.X) - Xm^2)
  X.cov = matrix(1/(p-1),p,p)
  diag(X.cov) = 1
  X.cov[1,-1] = X.cov[-1,1] = 0.5/sqrt(p-1)*(2-1/(p-1))
  
  X = (floor(log(1+exp(X.factor*Xraw)))-Xm)/X.sd
  S = (floor(log(1+exp(X.factor*Sraw)))-Xm)/X.sd
  
  return(list(X=X, S=S, Xraw = Xraw,X.cov = X.cov))
}

sim.cov.simple = function(N,p,q)
{
  X = matrix(rnorm(N*(p)),N,p)
  S = matrix(rnorm(N*q),N,q)
  
  return(list(X=X, S=S))
}

sim.SY.emr2 = function(Xraw, X.cov, n,m,M.coef,Y.coef, 
                        s.sd = 1, X.factor = 5, 
                        kterm = 100)
{
  N = nrow(Xraw)
  
  Smraw = (Xraw %*% M.coef)+matrix(rnorm(N*m),N,m)*s.sd
  Yseed = runif(N)
  
  Smraw.sd = sqrt(diag(t(M.coef)%*% X.cov %*% M.coef))
  prob.Sm = apply(pnorm(outer(log(exp(0:(kterm+1))-1)/X.factor,Smraw.sd,'/')),2,diff)
  Smm = apply((0:kterm)*prob.Sm,2,sum)
  Sm.sd = sqrt(apply((0:kterm)^2*prob.Sm,2,sum)-Smm^2)
  Sm = (floor(log(1+exp(X.factor*Smraw)))-rep(Smm,each=N))/rep(Sm.sd,each=N)
  Sm = Sm + 0.3
  Y = as.integer(Yseed <= expit(drop(as.matrix(Sm) %*% Y.coef)))
  
  return(list(S = Sm, Y = Y))
}

sim.SY.simple = function(X, n,m,M.coef,Y.coef, 
                       s.sd = 1)
{
  N = nrow(X)
  Sm = (X %*% M.coef)+matrix(rnorm(N*m),N,m)*s.sd
  Yseed = runif(N)
  Y = as.integer(Yseed <= expit(drop(as.matrix(Sm) %*% Y.coef)))
  
  return(list(S = Sm, Y = Y))
}

sim.auc = function(Y,X,S,M.coef, n, m)
{
  q = ncol(S)
  p = ncol(X)
  N = length(Y)
  M.cat = factor(M.coef)
  beta.grp = as.numeric(M.cat)
  ngrp = nlevels(M.cat)
  Xgrp = X %*% model.matrix(~-1 + M.cat)
  # 
  # sbeta = sum(M.coef== max(M.coef))
  # Xin = apply(X[,1:sbeta],1,sum)
  # Xout = apply(X[,-1:-sbeta], 1, sum)
  
  out= list()
  
  X.fit = glm(Y~Xgrp,family = binomial,subset = 1:n)
  W.fit = glm(Y~Xgrp+S[,1:m],family = binomial,subset = 1:n)
  out$beta.oracle = coef(X.fit)[c(1,beta.grp+1)]
  out$gamma.oracle = c(coef(W.fit)[c(1,beta.grp+1,ngrp+1+1:m)],rep(0,q-m))
  X.pred = expit(drop(out$beta.oracle[1] + X[-1:-n,] %*% out$beta.oracle[-1]))
  W.pred = expit(drop(out$gamma.oracle[1] + X[-1:-n,] %*% out$gamma.oracle[1+1:p]
                      + S[-1:-n, ] %*% out$gamma.oracle[-(1+0:p)]))
  null.pred = rep(mean(Y[1:n]),N-n)
  out$dev0 = logistic.dev(Y[-1:-n], null.pred)
  out$R2X = 1-logistic.dev(Y[-1:-n], X.pred)/logistic.dev(Y[-1:-n], null.pred)
  out$R2W = 1-logistic.dev(Y[-1:-n], W.pred)/logistic.dev(Y[-1:-n], null.pred)
  out$AUCX = auc(Y[-1:-n],X.pred)
  out$AUCW = auc(Y[-1:-n],W.pred)
  
  return(out)
}
