
rm(list = objects())

if(dir.exists("~/R-4.0.1/library"))
{
  .libPaths(c("~/R-4.0.1/library",.libPaths()))
  np = 20
}else{
  np = 8
}

pak.list = c("data.table", "bit64", "flexmix", "Matrix", "stats", "pROC",
             "splines", "randomForest")

for (pak in pak.list)
{
  yo = require(pak, character.only = T)
  if(!yo)
  {
    install.packages(pak,repos = "http://cran.us.r-project.org")
    require(pak, character.only = T)
  }
}

source("core/ate_np.R")
source("source/sim_gen_np.R")
source("core/cf_np.R")


# Find the parameters for different AUCs

S.gen.A = "beta.bias.gen"
S.gen.Y = "beta.bias.gen"
Y.gen = "binary.Y"
Y.par = list()

n = 500
N = 10000
Nrep = 1000
nfold = 5
nseq = 200

seq.ba = seq(1,7, length.out = nseq)
seq.by = seq(1,7, length.out = nseq)

grid.file = "simresult/np_beta_grid.rda"

if(file.exists(grid.file))
{
  load(grid.file)
}else{
  
  auc.ba = auc.by = rep(0,nseq)
  
  for (i in 1:nseq)
  {
    S.par.A = list(a0 = 1, a1 = seq.ba[i], b0 = seq.ba[i], b1 = 1)
    S.par.Y = list(a0 = 1, a1 = seq.by[i], b0 = seq.by[i], b1 = 1)
    
    dat = sim.gen.np(N, Y.gen = Y.gen, Y.par = Y.par, 
                     S.gen.A = S.gen.A, S.par.A = S.par.A, 
                     S.gen.Y = S.gen.Y, S.par.Y = S.par.Y)
    auc.by[i] = dat$imp.meas[2]
    auc.ba[i] = dat$imp.meas[1]
  }
  save(auc.ba, auc.by, file = grid.file)
}

ngrid = 5
auc.grid = seq(0.75, 0.95, length.out = ngrid)
auc.grid = c(0.80, 0.90, 0.95, 0.99, 0.999)
grid.ba = seq.ba[apply(abs(outer(auc.grid, auc.ba, '-')),1,which.min)]
grid.by = seq.by[apply(abs(outer(auc.grid, auc.by, '-')),1,which.min)]

library(doParallel)
library(doRNG)

cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)
if(dir.exists("~/R-4.0.1/library"))
{
  clusterEvalQ(cl, .libPaths(c("~/R-4.0.1/library",.libPaths())))
}

fit.fun = function(formula,...){coef(glm(formula,family = binomial,...))}

for(ia in 1:ngrid)
{  
  S.par.A = list(a0 = 1, a1 = grid.ba[ia], b0 = grid.ba[ia], b1 = 1)
  for (iy in 1:ngrid)
  {
    S.par.Y = list(a0 = 1, a1 = grid.by[iy], b0 = grid.by[iy], b1 = 1)
    
    setting = paste0("auc", round(auc.grid[ia]*100),
                     round(auc.grid[iy]*100))
    outfile = paste0("simresult/np_beta_bs_",
                      setting,"_v3.rda")
    if(file.exists(outfile))
      next
    
    registerDoRNG(seed = 531)
    sim.np <- foreach(irep = 1:Nrep,.packages = pak.list ) %dopar%
      {
        source("core/ate_np.R", local = T)
        source("source/sim_gen_np.R", local = T)
        source("core/cf_np.R", local = T)
        
        dat = sim.gen.np(N, Y.gen = Y.gen, Y.par = Y.par, 
                         S.gen.A = S.gen.A, S.par.A = S.par.A, 
                         S.gen.Y = S.gen.Y, S.par.Y = S.par.Y)
        W = cbind(dat$X, dat$S)
        
        foldid = c(rep_len(1:nfold, length.out = n)[sample(n)],
                   rep_len(1:nfold, length.out = N-n)[sample(N-n)])
        
        # Fit models over labeled data: splines
        #-----------------------------------------------------------
        
        # Shared models
        pi.bs = cf.bs(dat$A[1:n], dat$X, nfold, foldid,
                      fit.fun=fit.fun, link = expit, 
                      dev.fun = dev.binary.np)
        # pi.bs$meas
        mu1.bs = cf.bs(dat$Y[1:n], dat$X, nfold, foldid,
                       subset = dat$A[1:n]==1,
                       fit.fun=fit.fun, link = expit, 
                       dev.fun = dev.binary.np)
        # mu1.bs$meas
        mu0.bs = cf.bs(dat$Y[1:n], dat$X, nfold, foldid,
                       subset = dat$A[1:n]==0,
                       fit.fun=fit.fun, link = expit, 
                       dev.fun = dev.binary.np)
        # mu0.bs$meas
        
        # SSL models
        imp.A.bs = cf.bs(dat$A[1:n], W, nfold, foldid)
        # imp.A.bs$meas
        # imp.A1Y1.bs = cf.bs(dat$A[1:n]*dat$Y[1:n], W, nfold, foldid)
        imp.A1Y1.bs = cf.bs(dat$Y[1:n], W, nfold, foldid,
                            subset = dat$A[1:n]==1,
                            fit.fun=fit.fun, link = expit, 
                            dev.fun = dev.binary.np)
        # imp.A1Y1.bs$meas
        # imp.A0Y1.bs = cf.bs((1-dat$A[1:n])*dat$Y[1:n], W, nfold, foldid)
        imp.A0Y1.bs = cf.bs(dat$Y[1:n], W, nfold, foldid,
                            subset = dat$A[1:n]==0,
                            fit.fun=fit.fun, link = expit, 
                            dev.fun = dev.binary.np)
        # imp.A0Y1.bs$meas
        
        # Unsupervised
        #-----------------------------------------------------------
        SA = as.numeric(dat$S[,1]>=quantile(dat$S[,1], prob = 1-mean(dat$A[1:n])))
        SY = as.numeric(dat$S[,2]>=quantile(dat$S[,2], prob = 1-mean(dat$Y[1:n])))
        
        pi.bs = cf.bs(SA, dat$X, nfold, foldid,
                      fit.fun=fit.fun, link = expit,
                      dev.fun = dev.binary.np)
        mu1.bs = cf.bs(SY, dat$X, nfold, foldid,
                       subset = SA==1,
                       fit.fun=fit.fun, link = expit,
                       dev.fun = dev.binary.np)
        mu0.bs = cf.bs(SY, dat$X, nfold, foldid,
                       subset = SA==0,
                       fit.fun=fit.fun, link = expit,
                       dev.fun = dev.binary.np)
        
        UL.bs.fit = ate.SL(SY, SA,
                           mu1.bs$cf.pred, mu0.bs$cf.pred,
                           pi.bs$cf.pred, 1-pi.bs$cf.pred)
        
        # Estimators
        #-----------------------------------------------------------
        
        SL.bs.fit = ate.SL(dat$Y[1:n], dat$A[1:n], 
                           mu1.bs$cf.pred[1:n], mu0.bs$cf.pred[1:n], 
                           pi.bs$cf.pred[1:n], 1-pi.bs$cf.pred[1:n])
    
        SSL.bs.fit = ate.SSL(dat$Y[1:n], dat$A[1:n], 
                             mu1.bs$cf.pred, mu0.bs$cf.pred, 
                             pi.bs$cf.pred, 1-pi.bs$cf.pred, 
                             imp.A.bs$cf.pred, 
                             imp.A.bs$cf.pred*imp.A1Y1.bs$cf.pred,
                             (1-imp.A.bs$cf.pred)*imp.A0Y1.bs$cf.pred)
        
        
        true.ate = mean(dat$mu1-dat$mu0)
        
        # return the result
        return(list(true = true.ate, 
                    SSL.bs.est = SSL.bs.fit$est,
                    SSL.bs.se = SSL.bs.fit$se, 
                    SL.bs.est = SL.bs.fit$est,
                    SL.bs.se = SL.bs.fit$se 
                    ,UL.bs.est = UL.bs.fit$est,
                    UL.bs.se = UL.bs.fit$se
        ))
        
      }
    save(sim.np, file = outfile)
  }
}

stopCluster(cl)





