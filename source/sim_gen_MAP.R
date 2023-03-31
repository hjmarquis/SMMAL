require(data.table)
require(bit64)
require(flexmix)
require(Matrix)
require(stats)
require(pROC)

source("source/FUN_NLP_PheWAS_v2.R")

map.gen = function(Y, X, lam.note, a.cd, a.nlp, 
                   lam.cd0, lam.cd1, 
                   lam.nlp0, lam.npl1)
{
  # Testing parameters
  # Y = rbinom(100, 1, 0.5)
  # lam.cd0 = 0.5
  # lam.cd1 = 5
  # lam.nlp0 = 2
  # lam.npl1 = 10
  # lam.note = 50
  # a.cd = 1
  # a.nlp = 2
  
  # Function body
  
  n = length(Y)
  note = ceiling(rexp(n, 1/lam.note)+1)
  S.cd = rpois(n, a.cd*log(1+note)+Y*lam.cd1+(1-Y)*lam.cd0)
  S.nlp = rpois(n, a.nlp*log(1+note)+Y*lam.npl1+(1-Y)*lam.nlp0)
  
  MAPfit = MAP_PheWAS_main(dat.icd = list(ID=1:n, mat=data.frame(feature = S.cd)), 
                           dat.nlp = list(ID=1:n, mat=data.frame(feature = S.nlp)), 
                           dat.note = list(ID=1:n, mat=data.frame(feature = note)),
                           nm.phe = "feature",
                           p.icd = 0.001, n.icd = 10, 
                           p.nlp = 0.001, n.nlp = 10,
                           yes.con=FALSE, yes.nlp=FALSE)
  MAP.score = drop(MAPfit$res.all$feature$scores)
  
  # Test AUC
  # auc(Y, MAP.score, levels = c(0,1),
  #     direction = "<")
  
  return(list(S = MAP.score, 
              auc = auc(Y, MAP.score, levels = c(0,1),
                        direction = "<")))
}


pois.gen = function(Y, X, lam.0, lam.1)
{
  # Function body
  S = rpois(length(Y), Y*(lam.1-lam.0)+lam.0)

  
  return(list(S = S, 
              auc = auc(Y, S, levels = c(0,1),
                        direction = "<")))
}

beta.gen = function(Y, X, a0, b0, a1, b1)
{
  # Function body
  S = rbeta(length(Y), Y*(a1-a0)+a0, Y*(b1-b0)+b0)
  
  
  return(list(S = S, 
              auc = auc(Y, S, levels = c(0,1),
                        direction = "<")))
}

gaus.gen = function(Y, X, a, b, sd)
{
  # Function body
  S = rnorm(length(Y), a + b*Y, sd)
  
  return(list(S = S, Rsq = summary(lm(Y~S))$r.squared))
}


beta.bias.gen = function(Y, X, a0, b0, a1, b1)
{
  # Function body
  S = rbeta(length(Y), Y*(a1+X-a0)+a0,(1-Y)*(b0+X-b1)+b1)
  
  
  return(list(S = S, 
              auc = auc(Y, S, levels = c(0,1),
                        direction = "<")))
}


