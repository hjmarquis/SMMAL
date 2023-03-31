
source("source/sim_gen_MAP.R")

# Use the simulation from Rothe and Firpo 2018

RF18.pi = function(x, fmin= 0.4, fmax= 0.6)
{
  1-(1-fmax)/(1+((1-fmax)/(1-fmin) - 1)*x^2)
}

RF18.mu1 = function(x, fmin= 0.45, fmax= 0.65)
{
  1-(1-fmax)/(1+((1-fmax)/(1-fmin) - 1)*x^2)
}

RF18.mu0 = function(x, fmin= 0.35, fmax= 0.55)
{
  1-(1-fmax)/(1+((1-fmax)/(1-fmin) - 1)*(1-x)^2)
}

binary.Y = function(A, mu1, mu0)
{
  rbinom(length(A), 1, A*mu1+(1-A)*mu0)
}

gaus.Y = function(A, mu1, mu0, sd)
{
  rnorm(length(A), A*mu1+(1-A)*mu0, sd)
}

sim.gen.np = function(n, pi.fun = RF18.pi, 
                      mu1.fun = RF18.mu1, mu0.fun = RF18.mu0, 
                      Y.gen = "binary.Y", 
                      Y.par = c(),
                      S.gen.A = "map.gen",
                      S.par.A = list(lam.cd0 = 0.5, lam.cd1 = 5, 
                                   lam.nlp0 = 2, lam.npl1 = 10, 
                                   lam.note = 50, 
                                   a.cd = 1, a.nlp = 2), 
                      S.gen.Y = "map.gen", 
                      S.par.Y = list(lam.cd0 = 0.5, lam.cd1 = 5, 
                                     lam.nlp0 = 2, lam.npl1 = 10, 
                                     lam.note = 50, 
                                     a.cd = 1, a.nlp = 2))
{
  # Testing parameters
  # n = 100
  # map.par = list(lam.cd0 = 0.5, lam.cd1 = 5, 
  #                lam.nlp0 = 2, lam.npl1 = 10, 
  #                lam.note = 50, 
  #                a.cd = 1, a.nlp = 2)
  # pi.fun = RF18.pi
  # mu1.fun = RF18.mu1
  # mu0.fun = RF18.mu0
  
  # Body
  X = runif(n)
  
  pi = pi.fun(X)
  mu1 = mu1.fun(X)
  mu0 = mu0.fun(X)
  
  A = rbinom(n, 1, pi)
  Y = do.call(Y.gen, c(list(A=A, mu1 = mu1, mu0=mu0), Y.par))
  Sa = do.call(S.gen.A, c(list(Y=A,X=X), S.par.A))
  Sy = do.call(S.gen.Y, c(list(Y=Y,X=X), S.par.Y))
  
  return(list(Y=Y,A=A,X=matrix(X), S = cbind(Sa$S,Sy$S),
              pi = pi, mu1 = mu1, mu0 = mu0, 
              imp.meas = c(Sa[[2]], Sy[[2]])))
  
}
