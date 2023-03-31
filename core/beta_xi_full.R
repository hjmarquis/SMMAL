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
# xi = full.fit$xi[,1]
# beta = full.fit$beta[,1]
# lambda.list = full.fit$lambda
# 
# lambda = lambda * sqrt(0.95)
# 
# 
# xi = xi.lambda[,ilam-1]
# beta = beta.lambda[,ilam-1]
# init.xi1 = 0.1

beta.xi.full = function(X,S,trt,eta,n,
                        link, dlink, inv.link,
                        nlambda = 100, max.df = sqrt(n), 
                        lambda.list = NULL, 
                        maxit=1000, tol = 1e-7,
                        max.step = 1,
                        beta = rep(0,1+ncol(X)),
                        init.xi1 = 0.1,
                        xi = NULL)
{
  p = ncol(X)
  N = nrow(X)
  IX = cbind(1,X)
  W = cbind(1,X,S)
  pred.eta = expit(drop(W[-1:-n,]%*% eta))
  
  if(is.null(xi))
  {
    active = 1+which.max(abs(apply(W[1:n,-1], 2, cor,trt[1:n])))
    init.xi1 = init.xi1*sign(cor(W[1:n,active], trt[1:n]))
    xi0 = 0
    
    lp.beta = IX%*%beta
    pred.beta = link(lp.beta)
    dpred.beta = dlink(lp.beta[-1:-n])
    lp.xi = (W[,active]* init.xi1+xi0)
    pred.xi = link(lp.xi)
    dpred.xi = dlink(lp.xi[1:n])
    # initial beta
    for(iter in 1:(maxit*10))
    {
      score.xi = mean((pred.beta[1:n]*(trt[1:n]-pred.xi[1:n])))
      cont.xi= abs(score.xi)>tol
      if(cont.xi)
      {
        Hess = mean(dpred.xi*pred.beta[1:n])
        move = score.xi/Hess
        move.size = abs(move)
        if(move.size > max.step)
        {
          move = move/move.size*max.step
        }
        xi0 = xi0 + move
        lp.xi = lp.xi + move
        pred.xi = link(lp.xi)
        dpred.xi = dlink(lp.xi[1:n])
      }
      score = drop((pred.eta-pred.beta[-1:-n]*pred.xi[-1:-n]) %*% IX[-1:-n,])/N
      cont.beta = any(abs(score)>tol)
      if(cont.beta)
      {
        Hess = t(IX[-1:-n,])%*% (dpred.beta*pred.xi[-1:-n] * IX[-1:-n,])/N
        # Hess.imp = (((tmp.mf)^(-2) * t( %*% IX[1:n,])/n)  %*%
        #               (tmp.mf*((dlink(lp[1:n])*trt[1:n]) %*% IX[1:n,])/n
        #                - mean(link(lp[1:n])*trt[1:n])*(dlink(lp[1:n]) %*% IX[1:n,])/n)
        # )
        move = solve(Hess,score)
        move.size = sqrt(sum(move^2))
        if(move.size > max.step)
        {
          move = move/move.size*max.step
        }
        beta = beta + move
        lp.beta = IX%*%beta
        pred.beta = link(lp.beta)
        dpred.beta = dlink(lp.beta[-1:-n])
      }
      # print(xi0)
      # print(c(score.xi, max(abs(score))))
      if(!(cont.xi | cont.beta))
        break
    }
    
    if(iter == (maxit*10))
    {
      stop("Failed at initialization.")
    }

  
  
    # solution path 
    resn = (trt[1:n]-pred.xi[1:n])
    score.xi = c(drop((pred.beta[1:n]*resn) %*% W[1:n,1:(p+1)]),
      drop(resn %*% W[1:n,-1:-(p+1)]))/n
    
    lambda = max(abs(score.xi))
    if(is.null(lambda.list))
    {
      lambda.list = lambda* sqrt(0.95)^(1:nlambda - 1)
    }else{
      lambda.list = c(lambda, lambda.list)
      nlambda = length(lambda.list)
    }
    
    xi.lambda = sparseMatrix(1,1, x = xi0, 
                             dims = c(ncol(W),nlambda))
    beta.lambda = matrix(0, p+1, nlambda)
    active = c(1,active)
    xi = rep(0,ncol(W))
    xi[active] = c(xi0,init.xi1)
  }else
  {
    lp.beta = IX%*%beta
    pred.beta = link(lp.beta)
    dpred.beta = dlink(lp.beta[-1:-n])
    lp.xi = drop(W %*% xi)
    pred.xi = link(lp.xi)
    dpred.xi = dlink(lp.xi[1:n])
    active = c(1,1+which(xi[-1]!=0))
    if(is.null(lambda.list))
    {
      stop("Must provide lambda with initial xi.")
    }
    nlambda = length(lambda.list)
    xi.lambda = sparseMatrix(1,1, x = 0, 
                             dims = c(ncol(W),nlambda))
    beta.lambda = matrix(0, p+1, nlambda)
  }
  # load("other/test/test_dat.rda")
  # ilam = 80
  # xi = as.numeric(xi)
  
  for(ilam in 1:nlambda)
  {
    # print(ilam)
    lambda = lambda.list[ilam]
    for(iter in 1:maxit)
    {
      resn = (trt[1:n]-pred.xi[1:n])
      score.xi = c(drop((pred.beta[1:n]*resn) %*% W[1:n,1:(p+1)]),
                   drop(resn %*% W[1:n,-1:-(p+1)]))/n
      active = sort(unique(c(active,which(abs(score.xi) > lambda+tol))))
      score.xi.active = score.xi[active]
      score.xi.active[-1] = score.xi.active[-1] - sign(score.xi.active[-1])*lambda
      
      cont.xi = any(abs(score.xi.active)>tol)
      if(cont.xi)
      {
        activeX = active<= p+1
        Hess.active = rbind(t(W[1:n,active[activeX]])%*%
                              (pred.beta[1:n]*dpred.xi*W[1:n, active]),
                            t(W[1:n,active[!activeX]])%*%
                              (dpred.xi*W[1:n, active]))/n
        move <- tryCatch(solve(Hess.active,score.xi.active),
                         error = function(e) {return(NA)})
        if(any(is.na(move)))
        {
          iter = maxit
          break
        }
        move.size = sqrt(sum(move^2))
        if(move.size > max.step)
        {
          move = move/move.size*max.step
        }
        deactive = setdiff(which((sign(move)== sign(-xi[active])) & 
                                   (abs(move) >= abs(xi[active]))),
                           1)
        if(length(deactive)>0)
        {
          xi[active[deactive]] = 0
          move = move[-deactive]
          active = active[-deactive]
        }
        xi[active] = xi[active] + move
        lp.xi = lp.xi + drop(W[,active] %*% move)
        pred.xi = link(lp.xi)
        dpred.xi = dlink(lp.xi[1:n])
      }
      score =  drop((pred.eta-pred.beta[-1:-n]*pred.xi[-1:-n]) %*% IX[-1:-n,])/N
      cont.beta =  any(abs(score) > tol)
      if(cont.beta)
      {
        Hess = t(IX[-1:-n,])%*% (dpred.beta*pred.xi[-1:-n] * IX[-1:-n,])/N
        move <- tryCatch(solve(Hess,score),
                 error = function(e) {return(NA)})
        if(any(is.na(move)))
        {
          iter = maxit
          break
        }
        move.size = sqrt(sum(move^2))
        if(move.size > max.step)
        {
          move = move/move.size*max.step
        }
        beta =  beta + move
        lp.beta = IX%*%beta
        pred.beta = link(lp.beta)
        dpred.beta = dlink(lp.beta[-1:-n])
      }      
      # print(xi[active])
      # print(c(score.xi[1],
      # max(abs(score.xi))-lambda,max(abs(score))))
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

