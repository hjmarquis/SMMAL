cal.fit = function(trt,X,w, init,tol=1e-7,
                   maxit = 100, min.factor = 0.75,
                   ls.factor = 0.75, max.move = 1)
{
  n = nrow(X)
  X = cbind(1,X)
  oldscore = NULL
  
  for(k in 1:maxit) 
  {
    exp.lp = exp(-drop(X%*%init))
    bscore = apply((1-trt*(1+exp.lp))*w*X,2,mean)
    if(!is.null(oldscore))
      if(((sum(oldscore^2)*min.factor) <= sum(bscore^2)))
      {
        init = init+dinit
        dinit = dinit*ls.factor
        if(max(abs(dinit))<tol)
        {
          if(max(abs(oldscore)) > 1e-6)
            warning(paste("Algorithm stops in line-search. Target tol: ",
                          tol, ". Current tol: ", max(abs(oldscore)),
                          ". ", sep = ''))
          break
        }
        init = init - dinit
        next
      }
    oldscore = bscore
    bHess = t(X) %*% (trt*w*exp.lp*X) / n
    dinit = solve(bHess,bscore)
    if(all(abs(bscore)<tol))
      break
    # print(rbind(init,bscore,dinit))
    init = init - dinit
  }
  if(k >=maxit)
    stop("Numerical error when computing beta_delta")
  
  return(init)
}