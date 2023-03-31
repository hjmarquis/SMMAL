sim.Cox.logit = function(N,p,alpha,beta1,beta0,bash,tau,
                         zeta, s0rate, s1rate)
{
  X = matrix(rnorm(N*p),N,p)
  Xsq = scale(X[,1]^2)
  ps = 1/(1+exp(-drop(X %*% alpha[-1] + alpha[1]+Xsq*alpha.sq)))
  trt = rbinom(N,1,ps)
  
  pY1 = 1-exp(-tau*bash*exp(drop(X %*% beta1[-1]+beta1[1])))
  pY0 = 1-exp(-tau*bash*exp(drop(X %*% beta0[-1]+beta0[1])))
  pY = trt*pY1 + (1-trt)*pY0
  Y = rbinom(N,1,pY)
  
  rSt = trt*s1rate + (1-trt)*s0rate
  rSy = Y*s1rate + (1-Y)*s0rate
  S = cbind(rpois(N,rSt),
            rpois(N,rSy)) 
  
  impA = 1/(1+(1/ps-1)*dpois(S[,1],s0rate)/dpois(S[,1],s1rate))
  impY1 = 1/(1+(1/pY1-1)*dpois(S[,2],s0rate)/dpois(S[,2],s1rate))
  impY0 = 1/(1+(1/pY0-1)*dpois(S[,2],s0rate)/dpois(S[,2],s1rate))
  
  psR = 1/(1+exp(-drop(X %*% zeta[-1] + zeta[1])))
  R = rbinom(N,1,psR)
  Y[R==0] = NA
  trt[R==0] = NA
  
  return(list(data = data.frame(Y = Y, trt = trt, 
              X = X, Xsq = Xsq, S = S,
              R=R),
              ate = mean(pY1-pY0),
              or1 = pY1, or0 = pY0, 
              psA = ps, psR = psR, 
              impA = impA,
              impY1 = impY1,
              impY0 = impY0))
}


# ps.imp = 1/(1+(1/ps-1)*dpois(S[,1],s0rate)/dpois(S[,1],s1rate))
# apply((trt-ps.imp)*cbind(1,X,S),2,mean)
# boxplot(ps.imp~trt)
# auc(trt,ps.imp)
# 
# Y1.imp = 1/(1+(1/pY1-1)*dpois(S[,2],s0rate)/dpois(S[,2],s1rate))
# apply((Y[trt==1]-Y1.imp[trt==1])*cbind(1,X[trt==1,],S[trt==1,]),2,mean)
# apply((trt*Y-ps.imp*Y1.imp)*cbind(1,X,S),2,mean)
# boxplot(Y1.imp[trt==1]~Y[trt==1])
# auc(Y[trt==1],Y1.imp[trt==1])
# 
# mean(impA/ps*(impY1-pY1))
# mean((1-impA)/(1-ps)*(impY0-pY0))
# 
mean(R*(trt*Y-impA*impY1)/(ps*psR))
mean(R*pY1*(trt-impA)/(ps*psR))
# mean((trt*Y-impA*impY1)/(ps))
# mean(pY1*(trt-impA)/(ps))
# 
mean(R*((1-trt)*Y-(1-impA)*impY0)/((1-ps)*psR))
mean(R*pY0*((1-trt)-(1-impA))/((1-ps)*psR))
# mean(((1-trt)*Y-(1-impA)*impY0)/(1-ps))
# mean(pY0*((1-trt)-(1-impA))/(1-ps))


mean(R*(trt*Y-impAY1)/(ps*mean(R)))
mean(R*(trt*Y-impA*impY1)/(ps*mean(R)))
mean(R*pY1*(trt-impA)/(ps*mean(R)))
mean(R*((1-trt)*Y-impAY0)/((1-ps)*mean(R)))
mean(R*((1-trt)*Y-(1-impA)*impY0)/((1-ps)*mean(R)))
mean(R*pY0*((1-trt)-(1-impA))/((1-ps)*mean(R)))


mean(psR*((1-trt)*Y-(1-impA)*impY0)/((1-ps)*mean(R)))

apply(((1-trt)*Y-(1-impA)*impY0)*cbind(1,X),2,mean)
apply((trt*Y-impA*impY1)*cbind(1,X),2,mean)
