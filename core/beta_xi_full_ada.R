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
# max.df = sqrt(n)
# max.step = 1
# lambda.list = full.fit$lambda[-1]
# X = X[SSL.pos,]
# S = S[SSL.pos,]
# trt = trt[fold.pos]
# n = fold.nlab
# xi = full.fit$xi[,length(full.fit$lambda)]
# beta = full.fit$beta[,length(full.fit$lambda)]
# 
# lambda = lambda * sqrt(0.95)
# 
# 
# xi = xi.lambda[,ilam-1]
# beta = beta.lambda[,ilam-1]
# penalty.factor = rep(c(1,0,1),c(300,1,11))
xi = rep(0,ncol(W))
beta = rep(0,1+ncol(X))
penalty.factor = rep(1,ncol(X)+ncol(S))

beta.xi.full = function(X,S,trt,eta,n,
                        link, dlink, inv.link,
                        nlambda = 100, max.df = sqrt(n), 
                        lambda.list = NULL, 
                        maxit=1000, tol = 1e-7,
                        max.step = 1,
                        beta = rep(0,1+ncol(X)),
                        xi = rep(0,1+ncol(X)+ncol(S)), 
                        penalty.factor = rep(1,ncol(X)+ncol(S)))
{
  p = ncol(X)
  N = nrow(X)
  IX = cbind(1,X)
  W = cbind(1,X,S)
  pred.eta = expit(drop(W%*% eta))
  
  nlambda = ifelse(is.null(lambda.list), nlambda,
                   length(lambda.list)+1)
  
  xi.lambda = sparseMatrix(1,1,x=0,
                           dims = c(ncol(W),nlambda))
  beta.lambda = matrix(0, p+1, nlambda)
  # beta.lambda[,1] = beta
  active = c(1, 1+which(penalty.factor==0))
  # xi = c(inv.link(pred.xi),rep(0,ncol(W)-1))
  
  
  # load("other/test/test_dat.rda")
  # ilam = 80
  # xi = as.numeric(xi)
  
  for(ilam in 1:nlambda)
  {
    # print(ilam)
    if(ilam==1)
    {
      lambda = Inf
      lambda[penalty.factor==0] = 0 
    }else if(ilam == 2)
    {
      score.xi = drop((pred.lab*(trt[1:n]-link(lp.xi))) %*% W[1:n,])/n
      pos.pen = which(penalty.factor>0)
      lambda = max(abs(score.xi[1+pos.pen])/penalty.factor[pos.pen])
      if(is.null(lambda.list))
      {
        lambda.list = lambda* sqrt(0.95)^(1:nlambda - 1)
      }else{
        lambda.list = c(lambda, lambda.list)
        nlambda = length(lambda.list)
      }
      lambda = lambda.list[ilam]*penalty.factor
    }else{
      lambda = lambda.list[ilam]*penalty.factor
    }
    for(iter in 1:maxit)
    {
      lp = IX%*%beta
      pred.lab = link(lp[1:n])
      lp.xi = drop(W[1:n,]%*% xi)
      score.xi = drop((pred.lab*(trt[1:n]-link(lp.xi))) %*% W[1:n,])/n
      active = sort(unique(c(active,
                             1+which(abs(score.xi[-1]) > lambda+tol))))
      score.xi.active = score.xi[active]
      score.xi.active[-1] = (score.xi.active[-1] - 
                               sign(score.xi.active[-1])*lambda[active-1])
      
      cont.xi = any(abs(score.xi.active)>tol)
      if(cont.xi)
      {
        Hess.active = t(W[1:n,active])%*% (dlink(lp.xi)*pred.lab * W[1:n,active])/n
        move = solve(Hess.active,score.xi.active)
        deactive = setdiff(which((sign(move)== sign(-xi[active])) & 
                                   (abs(move) >= abs(xi[active]))),
                           1)
        move.size = sqrt(sum(move^2))
        if(move.size > max.step)
        {
          move = move/move.size*max.step
        }
        if(length(deactive)>0)
        {
          xi[active[deactive]] = 0
          move = move[-deactive]
          active = active[-deactive]
        }
        xi[active] = xi[active] + move
      }
      pred.xi = link(drop(W[-1:-n,]%*% xi))
      res = drop(pred.eta[-1:-n]-link(lp[-1:-n])*pred.xi)
      score = drop(res %*% IX[-1:-n,])/N
      cont.beta =  any(abs(score) > tol)
      if(cont.beta)
      {
        Hess = t(IX[-1:-n,])%*% (dlink(lp[-1:-n])*pred.xi * IX[-1:-n,])/N
        beta <-  tryCatch(beta + solve(Hess,score),
                          error = function(e) {return(NA)})
        if(any(is.na(beta)))
        {
          iter = maxit
          break
        }
      }      
      # print(c(score.xi[1],
      #         max(abs(score.xi[active[-1]])-lambda[active-1]),max(abs(score))))
      if(!(cont.xi | cont.beta))
        break
      
    }
    # print(c(ilam,score.xi[1],
    #         max(abs(score.xi))-lambda,max(abs(score))))
    if(iter == maxit)
    {
      warning(paste("Fail to converge at ",ilam,"th lambda.",
                    sep=''))
      return(list(lambda = lambda.list[1:(ilam-1)], 
                  beta = beta.lambda[,1:(ilam-1)], 
                  xi = xi.lambda[,1:(ilam-1)]))
    }
    active = which(xi!=0)
    xi.lambda[active,ilam] = xi[active]
    beta.lambda[,ilam] = beta
    
    if(sum(xi!=0) > max.df)
    {
      warning(paste("Reach max df ",ilam,"th lambda.",
                    sep=''))
      xi.lambda = xi.lambda[,1:ilam]
      beta.lambda = beta.lambda[,1:ilam]
      lambda.list = lambda.list[1:ilam]
      nlambda = ilam
      break
    }
  }
  # 
  # save(list=objects(), 
  #      file = "other/test/test_dat.rda")
  
  return(list(lambda = lambda.list, 
              beta = beta.lambda, 
              xi = xi.lambda))
}

