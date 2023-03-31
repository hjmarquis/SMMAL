
expit = function(x)
{
  1/(1+exp(-x))
}


dexpit = function(x)
{
  return(expit(x)*(1-expit(x)))
}


logit = function(x)
{
  log(x/(1-x))
}


logistic.dev = function(D,p,weights = rep(1,length(D)))
{
  if(sum(D==1)<=1)
    return(-(weights[D==1]*log(p[D==1,])+
               apply(weights[D==0]*log(1-p[D==0,]),2,sum))/sum(weights))
  if(sum(D==0)<=1)
    return(-(apply(weights[D==1]*log(p[D==1,]),2,sum)+
               weights[D==0]*log(1-p[D==0,]))/sum(weights))
  if(is.null(dim(p)))
    return(-(sum(weights[D==1]*log(p[D==1]))+
               sum(weights[D==0]*log(1-p[D==0])))/sum(weights))
  else
    return(-(apply(weights[D==1]*log(p[D==1,]),2,sum)+
               apply(weights[D==0]*log(1-p[D==0,]),2,sum))/sum(weights))
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