cf.beta = function(Y,X,S,gamma,nfold,Nfold, foldid)
{
  p = ncol(X)
  
  Ydag = expit(gamma$full[1]+
                 drop(X%*%gamma$full[1+1:p])
               + drop(S%*%gamma$full[-1:(-p-1)]))
  full = coef(glm(Ydag~X,family = binomial))
  full[is.na(full)] = 0

  if(nfold<=2)
  {
    return(full)
  }
  
  cf = vector("list",nfold)
  for(ifold in 1:nfold)
  {
    Ydag = expit(gamma$cf[1,ifold]+
                   drop(X%*%gamma$cf[1+1:p,ifold])
                 + drop(S%*%gamma$cf[-1:(-p-1),ifold]))
    cf[[ifold]] = matrix(0,length(full),nfold)
    for(jfold in 1:Nfold)
    {
      tmpcoef = coef(glm(Ydag~X,family = binomial,
                                     subset = foldid!=jfold))
      tmpcoef[is.na(tmpcoef)] = 0
      cf[[ifold]][,jfold] = tmpcoef
    }
  }
  
  return(list(full=full,cf=cf))
}