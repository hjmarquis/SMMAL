
expit = function(x)
{
  1/(1+exp(-x))
}


logit = function(x)
{
  log(x/(1-x))
}


logistic.dev = function(Y,lp,weights = rep(1,length(Y)))
{
  if(is.null(dim(lp)))
    return(sum(weights*(log(1+exp(lp)) - Y*lp))/sum(weights))
  else
    return(apply(weights*(log(1+exp(lp)) - Y*lp),2,sum)/sum(weights))
}

logistic.R2 = function(D,lp)
{
  dev = logistic.dev(D,expit(lp))
  dev0 = logistic.dev(D,rep(mean(D),length(D)))
  return(1-dev/dev0)
}

get.var = function(varmat,testx)
{
  apply(testx,1,function(x,V) drop(t(x)%*% V %*%x), V=varmat)
}