semi.eff.ate = function(R,Y,trt, or1, or0, 
                        psA1,psA0, psR, 
                        impA1, impA0, impY1,impY0)
{
  infl = or1 + impA1/psA1*(impY1-or1) - or0 - impA0/psA0*(impY0-or0)
  pos.label = which(R==1)
  infl[pos.label] = (infl[pos.label]+
    (trt[pos.label]*Y[pos.label]-impA1[pos.label]*impY1[pos.label]
     -or1[pos.label]*(trt[pos.label]-impA1[pos.label]))/
      (psA1[pos.label]*psR[pos.label])
    -((1-trt[pos.label])*Y[pos.label]-impA0[pos.label]*impY0[pos.label]
      -or0[pos.label]*((1-trt[pos.label])-impA0[pos.label]))/
      (psA0[pos.label]*psR[pos.label]))

  return(list(ate = mean(infl),
              se = sqrt(var(infl)/length(R))))
}

super.ate = function(R,Y,trt, or1, or0, 
                        psA1,psA0,psR)
{
  infl = or1 - or0
  pos.label = which(R==1)
  infl[pos.label] = (infl[pos.label] + trt[pos.label]/psA1[pos.label]/psR[pos.label]*(Y[pos.label]-or1[pos.label])
           - (1-trt[pos.label])/psA0[pos.label]/psR[pos.label]*(Y[pos.label]-or0[pos.label]))
 
  return(list(ate = mean(infl),
              se = sqrt(var(infl)/length(R))))
}