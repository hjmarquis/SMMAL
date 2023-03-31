expit = function(x)
{
  return(1/(1+exp(-x)))
}

dev.binary = function(y, pi)
{
  if(!is.null(dim(pi)))
  {
    return(-apply(y*log(pi)+(1-y)*log(1-pi),2,mean))
  }else{
    return(-mean(y*log(pi)+(1-y)*log(1-pi)))
  }
}

dev.binary.np = function(y, pi)
{
  if(!is.null(dim(pi)))
  {
    return(apply(pi,2,dev.binary.np.uni,
                  y=y))
  }else{
    return(dev.binary.np.uni(y,pi))
  }
}

dev.binary.np.uni = function(y,pi)
{
  pi = pmin(1,pmax(0, pi))
  if(any(y==1 & pi==0) | any(y==0 & pi==1))
  {
    return(Inf)
  }
  loglik = rep(NA, length(y))
  loglik[y==1] = log(pi[y==1])
  loglik[y==0] = log(1-pi[y==0])
  return(-mean(loglik))
}

dev.ls = function(y, mu, wgt = rep(1,length(y)))
{
  if(!is.null(dim(pi)))
  {
    return(apply((y - mu)^2*wgt,2,sum)/sum(wgt))
  }else{
    return(sum((y - mu)^2*wgt)/sum(wgt))
  }
}

dev.logistic = function(y, lp, wgt = rep(1,length(y)))
{
  if(!is.null(dim(lp)))
  {
    return(apply((log(1+exp(-lp))+(1-y)*lp)*wgt,2,sum)/sum(wgt))
  }else{
    return(sum((log(1+exp(-lp))+(1-y)*lp)*wgt)/sum(wgt))
  }
}


dev.rcal = function(y, lp, wgt = rep(1,length(y)))
{
  if(!is.null(dim(lp)))
  {
    return(apply(((1-y)*lp+y*exp(-lp))*wgt,2,sum) / sum(wgt))
  }else{
    return(sum(((1-y)*lp+y*exp(-lp))*wgt)/sum(wgt))
  }
}