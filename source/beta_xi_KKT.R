# 
# load("other/test/test_dat.rda")
# 
# # 
# # rm(list = setdiff(objects(), 
# #                   c("X","S","trt","eta","n")))
# # source("source/utility.R")
# maxit = 1000
# tol = 1e-7
# link = expit
# dlink = dexpit
# inv.link = logit
# ilam = 10
# eta = eta$full
# 
# beta = beta.xi.1$beta[,ilam]
# xi = beta.xi.1$xi[,ilam]
# lambda = beta.xi.1$lambda[ilam]


beta.xi.KKT= function(X,S,trt,eta,n,
                        link, dlink, inv.link,
                        beta, xi, lambda,  
                        maxit=1000, tol = 1e-7)
{
  # Check feasibility
  p = ncol(X)
  N = nrow(X)
  IX = cbind(1,X)
  W = cbind(1,X,S)
  pred.eta = expit(drop(W%*% eta))
  
  lp = IX%*%beta
  pred.lab = link(lp[1:n])
  pred = link(lp[-1:-n])
  lp.xi = drop(W%*% xi)
  res.xi = (trt[1:n]-link(lp.xi[1:n]))
  score.xi = drop((pred.lab*res.xi) %*% W[1:n,])/n
  
  pred.xi = link(drop(W[-1:-n,]%*% xi))
  res = drop(pred.eta[-1:-n]-link(lp[-1:-n])*pred.xi)
  score = drop(res %*% IX[-1:-n,])/N
  
  if(any(abs(score.xi) > lambda+tol) | any(abs(score) > tol)
     | (abs(score.xi[1]) > tol))
    stop("Solution violates the feasibility.")
  
  active.xi = setdiff(which(xi!=0), 1)
  active.score = which(abs(score.xi) > lambda-tol)
  
  # print(active.xi)
  # print(active.score)
  
  if((length(active.xi)) ==0)
    return(NULL)
  
  if((length(active.xi)>length(active.score)))
    stop("KKT condition not satisfied")
  
  partial.xi = -(t(W[1:n,active.score])%*% 
    (dlink(lp.xi[1:n])*pred.lab * W[1:n,active.xi])/n)
  partial.xi0 = -(t(W[1:n,active.score])%*% 
                    (dlink(lp.xi[1:n])*pred.lab)/n)
  partial.beta = (t(W[1:n,active.score])%*% 
                    (dlink(lp[1:n])*res.xi * IX[1:n,])/n)
  
  A1 = (t(IX[-1:-n,])%*% 
          (dlink(lp.xi[-1:-n])*pred * 1)/N)
  B1 = (t(IX[-1:-n,])%*% 
          (dlink(lp[-1:-n])*pred.xi * IX[-1:-n,])/N)
  C1 = - (t(IX[-1:-n,])%*% 
            (dlink(lp.xi[-1:-n])*pred * W[-1:-n,active.xi])/N)
  
  A2 = -mean(dlink(lp.xi[1:n])*pred.lab)
  B2 = apply(dlink(lp[1:n])*res.xi * IX[1:n,],2,mean)
  C2 = ((dlink(lp.xi[1:n])*pred.lab) %*% W[1:n,active.xi])/n
  
  Hess = (partial.xi + cbind(partial.xi0,partial.beta) %*% 
            solve(rbind(cbind(A1,B1),c(A2,B2)),rbind(C1,C2)))
  proj = solve(Hess, -sign(xi[active.xi]))
  if(any((proj * sign(score.xi[active.score]))<0))
    stop("KKT condition not satisfied")
  
  return(NULL)
}