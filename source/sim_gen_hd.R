
source("source/sim_gen_MAP.R")
source("core/utility.R")

# Use the simulation from Rothe and Firpo 2018

binary.Y = function(A, mu1, mu0)
{
  rbinom(length(A), 1, A*mu1+(1-A)*mu0)
}

X.iid = function(n, p)
{
  matrix(rnorm(n*p), n,p)
}

X.exch = function(n,p,rho)
{
  matrix(rnorm(n*p), n,p)*sqrt(1-rho) + rnorm(n)*sqrt(rho)
}

X.ar = function(n,p,coef)
{
  x = matrix(rnorm(n*p), n,p)
  for (j in 2:p) {
    x[,j] = drop(x[,j:(j+1-length(coef))] %*% coef)
  }
  return(x)
}

lp.correct = function(X,coef)
{
  return(drop(X %*% coef[-1])+coef[1])
}

lp.wrong = function(X, coef1, coef2)
{
  lp1 = drop(X %*% coef1[-1])
  lp2 = drop(X %*% coef2)
  return(lp1+coef1[1] + lp1*lp2)
}

sim.gen.hd = function(n, p, 
                      X.gen = "X.iid", par.X = list(),
                      lp.A = "lp.correct", par.lp.A, 
                      lp.Y = "lp.correct", par.lp.Y1, par.lp.Y0, 
                      link.A = "expit", par.link.A = list(), 
                      link.Y = "expit", par.link.Y = list(), 
                      Y.gen = "binary.Y", Y.par = list(), 
                      S.gen.A = "map.gen",
                      S.par.A = list(lam.cd0 = 1, lam.cd1 = 5, 
                                   lam.nlp0 = 4, lam.npl1 = 10, 
                                   lam.note = 50, 
                                   a.cd = 1, a.nlp = 2), 
                      S.gen.Y = "map.gen", 
                      S.par.Y = list(lam.cd0 = 1, lam.cd1 = 5, 
                                     lam.nlp0 = 4, lam.npl1 = 10, 
                                     lam.note = 50, 
                                     a.cd = 1, a.nlp = 2))
{
  # Testing parameters
  # n = 100
  # p=50
  # coef = c(0,1,rep(0,p-1))
  # X.gen = "X.iid"
  # par.X = list()
  # lp.A = "lp.correct"
  # lp.Y = "lp.correct"
  # par.lp.A = par.lp.Y1 = par.lp.Y0 = list(coef = coef) 
  # link.A = "expit"
  # par.link.A = list()
  # link.Y = "expit"
  # par.link.Y = list()
  # Y.gen = "binary.Y"
  # Y.par = list()
  # S.gen.A = "map.gen"
  # S.par.A = list(lam.cd0 = 0.5, lam.cd1 = 5, 
  #                lam.nlp0 = 2, lam.npl1 = 10, 
  #                lam.note = 50, 
  #                a.cd = 1, a.nlp = 2) 
  # S.gen.Y = "map.gen"
  # S.par.Y = list(lam.cd0 = 0.5, lam.cd1 = 5, 
  #                lam.nlp0 = 2, lam.npl1 = 10, 
  #                lam.note = 50, 
  #                a.cd = 1, a.nlp = 2)
  
  # Body
  X = do.call(X.gen, c(list(n=n,p=p), par.X))
  
  pi.lp = do.call(lp.A, c(list(X=X), par.lp.A))
  pi = do.call(link.A, c(list(x = pi.lp), par.link.A))
  
  
  mu1.lp = do.call(lp.Y, c(list(X=X), par.lp.Y1))
  mu0.lp = do.call(lp.Y, c(list(X=X), par.lp.Y0))
  mu1 = do.call(link.Y, c(list(x = mu1.lp), par.link.Y))
  mu0 = do.call(link.Y, c(list(x = mu0.lp), par.link.Y))
  
  A = rbinom(n, 1, pi)
  Y = do.call(Y.gen, c(list(A=A, mu1 = mu1, mu0=mu0), Y.par))
  Sa = do.call(S.gen.A, c(list(Y=A, X=pnorm(X[,1,drop = FALSE])), S.par.A))
  Sy = do.call(S.gen.Y, c(list(Y=Y, X=pnorm(X[,1,drop = FALSE])), S.par.Y))
  
  return(list(Y=Y,A=A,X=X, S = cbind(Sa$S,Sy$S),
              pi = pi, mu1 = mu1, mu0 = mu0, 
              imp.meas = c(Sa[[2]], Sy[[2]])))
  
}
