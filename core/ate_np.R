
ate.SL = function(Y, A, mu1, mu0, pi1,pi0, min.pi = 0.05, max.pi = 0.95)
{
  pi1 = pmin(max.pi, pmax(min.pi, pi1))
  pi1 = pi1 * mean(A)/mean(pi1)
  pi0 = pmin(max.pi, pmax(min.pi, pi0))
  pi0 = pi0 * mean(A)/mean(pi0)
  
  infl = mu1 + A/pi1*(Y-mu1) - mu0 - (1-A)/pi0*(Y-mu0)
  return(list(est = mean(infl),
              se = sd(infl)/sqrt(length(Y))))
}

ate.SSL = function(Y,A, mu1, mu0, pi1,pi0, imp.A,
                   imp.A1Y1, imp.A0Y1, min.pi = 0.05, max.pi = 0.95)
{
  n = length(Y)
  N = length(pi1)
  rho.inv = N/n
  
  pi1 = pmin(max.pi, pmax(min.pi, pi1))
  pi1 = pi1 * mean(A)/mean(pi1)
  pi0 = pmin(max.pi, pmax(min.pi, pi0))
  pi0 = pi0 * mean(A)/mean(pi0)
  
  infl = (
    (mu1 + imp.A1Y1/pi1 - imp.A*mu1/pi1) - 
    (mu0 + imp.A0Y1/pi0 - (1-imp.A)*mu0/pi0)
    )
  infl[1:n] = infl[1:n] + rho.inv*(
    (A*Y - imp.A1Y1[1:n])/pi1[1:n]
    -((1-A)*Y - imp.A0Y1[1:n])/pi0[1:n]
    - (A - imp.A[1:n])*(mu1[1:n]/pi1[1:n] + mu0[1:n]/pi0[1:n])
  )
  
  return(list(est = mean(infl),
              se = sd(infl)/sqrt(N)))
}