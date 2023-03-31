

sim.data.iid = function(N,p,q,m,s,M.coef,Y.coef, 
                        m.sd = sqrt(0.5), s.mu = 0.5, s.sd = 1,
                        s.hid = 0)
{
  X = matrix(rnorm(N*p),N,p)
  S = matrix(rnorm(N*q),N,q)
  Yseed = runif(N)
  
  MX = (X %*% M.coef)
  Sm = MX + S[,1:m]*m.sd 
  Y = as.integer(Yseed <= pnorm(drop(Sm %*% Y.coef+rnorm(N,sd=s.hid))))
  Ss = matrix(S[,m+1:s]*s.sd,N,s) + Y*s.mu
  S[,1:m] = Sm
  S[,m+1:s] = Ss
  
  # only works for symmetric M.coef
  sbeta = sum(M.coef!=0)
  Xin = apply(X[,1:sbeta],1,sum)
  Xout = apply(X[,-1:-sbeta], 1, sum)
  
  X.fit = glm(Y~Xin+Xout,family = binomial,subset = 1:n)
  W.fit = glm(Y~Xin+Xout+Sm + Ss,family = binomial,subset = 1:n)
  beta.oracle = rep(coef(X.fit),c(1,sbeta,p-sbeta))
  gamma.oracle = c(rep(coef(W.fit),c(1,sbeta,p-sbeta,rep(1,s+m))),rep(0,q-s-m))
  X.pred = expit(drop(beta.oracle[1] + X[-1:-n,] %*% beta.oracle[-1]))
  W.pred = expit(drop(gamma.oracle[1] + X[-1:-n,] %*% gamma.oracle[1+1:p]
                + S[-1:-n, ] %*% gamma.oracle[-(1+0:p)]))
  null.pred = rep(mean(Y[1:n]),N-n)
  R2X = 1-logistic.dev(Y[-1:-n], X.pred)/logistic.dev(Y[-1:-n], null.pred)
  R2W = 1-logistic.dev(Y[-1:-n], W.pred)/logistic.dev(Y[-1:-n], null.pred)
  AUCX = auc(Y[-1:-n],X.pred)
  AUCW = auc(Y[-1:-n],W.pred)
  
  return(list(X=X,S=S,Y=Y,
              beta.oracle=beta.oracle,
              gamma.oracle = gamma.oracle,
              R2X = R2X, R2W = R2W,
              AUCX = AUCX, AUCW = AUCW))
}

sim.data.exch = function(N,p,q,m,s,M.coef,Y.coef, 
                         m.sd = sqrt(0.5), s.mu = 0.5, s.sd = 1,
                         s.hid = 0)
{
  U = rnorm(N)/sqrt(p)
  X = matrix(rnorm(N*p),N,p)*sqrt(1-1/p)+U 
  S = matrix(rnorm(N*q),N,q)
  Yseed = runif(N)
  
  MX = (X %*% M.coef)
  Sm = MX + S[,1:m]*m.sd 
  Y = as.integer(Yseed <= pnorm(drop(Sm %*% Y.coef)+rnorm(N,sd=s.hid)))
  Ss = matrix(S[,m+1:s]*s.sd,N,s) + Y*s.mu
  S[,1:m] = Sm
  S[,m+1:s] = Ss
  
  # only works for symmetric M.coef
  sbeta = sum(M.coef!=0)
  Xin = apply(X[,1:sbeta],1,sum)
  Xout = apply(X[,-1:-sbeta], 1, sum)
  
  X.fit = glm(Y~Xin+Xout,family = binomial,subset = 1:n)
  W.fit = glm(Y~Xin+Xout+Sm + Ss,family = binomial,subset = 1:n)
  beta.oracle = rep(coef(X.fit),c(1,sbeta,p-sbeta))
  gamma.oracle = c(rep(coef(W.fit),c(1,sbeta,p-sbeta,rep(1,s+m))),rep(0,q-s-m))
  X.pred = expit(drop(beta.oracle[1] + X[-1:-n,] %*% beta.oracle[-1]))
  W.pred = expit(drop(gamma.oracle[1] + X[-1:-n,] %*% gamma.oracle[1+1:p]
                      + S[-1:-n, ] %*% gamma.oracle[-(1+0:p)]))
  null.pred = rep(mean(Y[1:n]),N-n)
  R2X = 1-logistic.dev(Y[-1:-n], X.pred)/logistic.dev(Y[-1:-n], null.pred)
  R2W = 1-logistic.dev(Y[-1:-n], W.pred)/logistic.dev(Y[-1:-n], null.pred)
  AUCX = auc(Y[-1:-n],X.pred)
  AUCW = auc(Y[-1:-n],W.pred)
  
  return(list(X=X,S=S,Y=Y,
              beta.oracle=beta.oracle,
              gamma.oracle = gamma.oracle,
              R2X = R2X, R2W = R2W,
              AUCX = AUCX, AUCW = AUCW))
}


sim.data.emr = function(N,p,q,m,s,M.coef,Y.coef, 
                        m.sd = sqrt(0.5), s.mu = 0.5, s.sd = 1,
                        X.factor = 3,
                        s.hid = 0)
{
  X.sd = sqrt(2/(1-1/base::pi)) * X.factor
  X.center = X.sd/sqrt(2*base::pi)
  
  U = rnorm(N)/sqrt(p)
  Xraw = matrix(rnorm(N*p),N,p)*sqrt(1-1/p)+U 
  Sraw = matrix(rnorm(N*q),N,q)
  Yseed = runif(N)
  
  X = round(apply(Xraw,2,pmax,0)*X.sd-X.center)/X.factor
  S = round(apply(Sraw,2,pmax,0)*X.sd-X.center)/X.factor
  MX = (X %*% M.coef)*X.factor
  Sm =  round(MX +apply(matrix(Sraw[,1:m],N),2,pmax,0)*X.sd*m.sd-X.center*m.sd)/X.factor
  Y = as.integer(Yseed <= pnorm(drop(Sm %*% Y.coef)+rnorm(N,sd=s.hid)))
  Ss = round(apply(matrix(Sraw[,m+1:s]*s.sd*X.sd,N,s) + Y*(s.mu*X.factor+X.center),2,pmax,0))/X.factor
  S[,1:m] = Sm
  S[,m+1:s] = Ss
  
  # only works for symmetric M.coef
  sbeta = sum(M.coef== max(M.coef))
  Xin = apply(X[,1:sbeta],1,sum)
  Xout = apply(X[,-1:-sbeta], 1, sum)
  
  X.fit = glm(Y~Xin+Xout,family = binomial,subset = 1:n)
  W.fit = glm(Y~Xin+Xout+Sm + Ss,family = binomial,subset = 1:n)
  beta.oracle = rep(coef(X.fit),c(1,sbeta,p-sbeta))
  gamma.oracle = c(rep(coef(W.fit),c(1,sbeta,p-sbeta,rep(1,s+m))),rep(0,q-s-m))
  X.pred = expit(drop(beta.oracle[1] + X[-1:-n,] %*% beta.oracle[-1]))
  W.pred = expit(drop(gamma.oracle[1] + X[-1:-n,] %*% gamma.oracle[1+1:p]
                      + S[-1:-n, ] %*% gamma.oracle[-(1+0:p)]))
  null.pred = rep(mean(Y[1:n]),N-n)
  R2X = 1-logistic.dev(Y[-1:-n], X.pred)/logistic.dev(Y[-1:-n], null.pred)
  R2W = 1-logistic.dev(Y[-1:-n], W.pred)/logistic.dev(Y[-1:-n], null.pred)
  AUCX = auc(Y[-1:-n],X.pred)
  AUCW = auc(Y[-1:-n],W.pred)
  
  return(list(X=X,S=S,Y=Y,
              beta.oracle=beta.oracle,
              gamma.oracle = gamma.oracle,
              R2X = R2X, R2W = R2W,
              AUCX = AUCX, AUCW = AUCW))
}


sim.data.emr2 = function(N,n,p,q,m,M.coef,Y.coef, 
                        s.sd = 1, X.factor = 5, 
                        kterm = 100)
{
  U = rnorm(N)/sqrt(p-1)
  Xraw = matrix(rnorm(N*(p-1)),N,p-1)*sqrt(1-1/(p-1))+U 
  X1raw = Xraw %*% rep(0.5/sqrt(p-1),p-1) + rnorm(N)*sqrt(.5)
  Xraw = cbind(X1raw,Xraw)
  Sraw = matrix(rnorm(N*(q-1)),N,q-1)
  Smraw = (Xraw %*% M.coef)+matrix(rnorm(N*m),N,m)*s.sd
  Yseed = runif(N)
  
  prob.X = diff(pnorm(log(exp(0:(kterm+1))-1)/X.factor))
  Xm = sum((0:kterm)*prob.X)
  X.sd = sqrt(sum((0:kterm)^2*prob.X) - Xm^2)
  X.cov = matrix(1/(p-1),p,p)
  diag(X.cov) = 1
  X.cov[1,-1] = X.cov[-1,1] = 0.5/sqrt(p-1)*(2-1/(p-1))
  Smraw.sd = sqrt(diag(t(M.coef)%*% X.cov %*% M.coef))
  prob.Sm = apply(pnorm(outer(log(exp(0:(kterm+1))-1)/X.factor,Smraw.sd,'/')),2,diff)
  Smm = apply((0:kterm)*prob.Sm,2,sum)
  Sm.sd = sqrt(apply((0:kterm)^2*prob.Sm,2,sum)-Smm^2)
  
  X = (floor(log(1+exp(X.factor*Xraw)))-Xm)/X.sd
  S = cbind(
    (floor(log(1+exp(X.factor*Smraw)))-rep(Smm,each=N))/rep(Sm.sd,each=N),
    (floor(log(1+exp(X.factor*Sraw)))-Xm)/X.sd)
  Y = as.integer(Yseed <= expit(drop(as.matrix(S[,1:m]) %*% Y.coef)))
    
  out = list(X=X,S=S,Y=Y)
  if(n>0)
  {
    # only works for symmetric M.coef
    M.cat = factor(M.coef)
    beta.grp = as.numeric(M.cat)
    ngrp = nlevels(M.cat)
    Xgrp = X %*% model.matrix(~-1 + M.cat)
    # 
    # sbeta = sum(M.coef== max(M.coef))
    # Xin = apply(X[,1:sbeta],1,sum)
    # Xout = apply(X[,-1:-sbeta], 1, sum)

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
  }
  
  return(out)
}