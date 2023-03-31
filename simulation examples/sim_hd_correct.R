
rm(list = objects())

source("core/ate_np.R")
source("source/sim_gen_hd.R")
source("core/cf_hd.R")

S.gen.A = "beta.bias.gen"
S.par.A = list(a0 = 1, a1 = 2, b0 = 2, b1 = 1)

S.gen.Y = "beta.bias.gen"
S.par.Y = list(a0 = 1, a1 = 2, b0 = 2, b1 = 1)

if(dir.exists("~/R-4.0.1/library"))
{
  .libPaths(c("~/R-4.0.1/library",.libPaths()))
  np = 20
}else{
  np = 8
}

pak.list = c("data.table", "bit64", "flexmix", "Matrix", "stats", "pROC",
                 "glmnet", "RCAL", "doParallel", "doRNG")

for (pak in pak.list)
{
  yo = require(pak, character.only = T)
  if(!yo)
  {
    install.packages(pak,repos = "http://cran.us.r-project.org")
    require(pak, character.only = T)
  }
}

cl = makeCluster(getOption("cl.cores", np))
registerDoParallel(cl)
registerDoRNG(seed = 531)
if(dir.exists("~/R-4.0.1/library"))
{
  clusterEvalQ(cl, .libPaths(c("~/R-4.0.1/library",.libPaths())))
}

n = 500
N = 10000
p = 500
coef.a = c(0,0.5, 0.25,0.125,0, rep(0,p-4))
coef.b1 = c(0.2,0.5, 0.25,0.125,0, rep(0,p-4))/2
coef.b0 = -c(0.2,0.5, 0.25,0.125,0, rep(0,p-4))/2
Nrep = 500
nfold = 5
nfold.ds = 3

sim.hd <- foreach(irep = 1:Nrep,.packages = pak.list ) %dopar%
  {
    source("core/ate_np.R", local = T)
    source("source/sim_gen_hd.R", local = T)
    source("core/cf_hd.R", local = T)
    source("core/ds_hd.R", local = T)
    
    dat = sim.gen.hd(N, p, 
                     par.lp.A = list(coef=coef.a),
                     par.lp.Y1 = list(coef=coef.b1),
                     par.lp.Y0 = list(coef=coef.b0),
                     S.gen.A = S.gen.A, S.par.A = S.par.A, 
                     S.gen.Y = S.gen.Y, S.par.Y = S.par.Y)
    # mean(dat$Y[dat$A==1]) - mean(dat$Y[dat$A==0])
    # mean(dat$mu1) - mean(dat$mu0)
    # summary(dat$pi)
    rnd.n = sample(n)
    rnd.N = sample(N-n)
    
    
    W = cbind(dat$X, dat$S)
    
    dsid = c(rep_len(1:2, length.out = n)[rnd.n],
             rep_len(1:2, length.out = N-n)[rnd.N])
    
    foldid = c(rep_len(1:nfold, length.out = n)[rnd.n],
               rep_len(1:nfold, length.out = N-n)[rnd.N])
    
    save(dat, rnd.n, rnd.N, file = paste0("tmpdat/hd_correct_dat_rep",
                                          irep, ".rda"))
    
    # SMMAL
    #-----------------------------------------------------------
    
    # Shared models
    pi.lasso = cf.lasso(dat$A[1:n], dat$X, nfold, foldid,
                        dcf = TRUE)
    
    # pi.bs$meas
    mu1.lasso = cf.lasso(dat$Y[1:n], dat$X, nfold, foldid,
                         subset = dat$A[1:n]==1,
                         dcf = TRUE)
    # mu1.bs$meas
    mu0.lasso = cf.lasso(dat$Y[1:n], dat$X, nfold, foldid,
                         subset = dat$A[1:n]==0,
                         dcf = TRUE)
    
    # mu0.bs$meas
    pi1.rcal = cf.rcal(dat$A[1:n], dat$X, nfold, foldid,
                       weight = mu1.lasso$cf.pred *(1-mu1.lasso$cf.pred),
                       dcf.weight = mu1.lasso$dcf.pred*(1-mu1.lasso$dcf.pred))
    pi0.rcal = cf.rcal(1-dat$A[1:n], dat$X, nfold, foldid,
                       weight = mu0.lasso$cf.pred *(1-mu0.lasso$cf.pred),
                       dcf.weight = mu0.lasso$dcf.pred*(1-mu0.lasso$dcf.pred))
    mu1.rw = cf.lasso(dat$Y[1:n], dat$X, nfold, foldid,
                      subset = dat$A[1:n]==1,
                      weight = (1/pi.lasso$cf.pred -1),
                      dcf.weight = (1/pi.lasso$dcf.pred -1))
    mu0.rw = cf.lasso(dat$Y[1:n], dat$X, nfold, foldid,
                      subset = dat$A[1:n]==0,
                      weight = (1/(1-pi.lasso$cf.pred) -1),
                      dcf.weight = (1/(1-pi.lasso$dcf.pred) -1))
    
    # SSL models
    imp.A.lasso = cf.lasso(dat$A[1:n], W, nfold, foldid)
    
    imp.A1Y1.lasso = cf.lasso(dat$Y[1:n], W, nfold, foldid,
                              subset = dat$A[1:n]==1)
    
    imp.A0Y1.lasso = cf.lasso(dat$Y[1:n], W, nfold, foldid,
                              subset = dat$A[1:n]==0)
    
    
    
    # Estimators
    #-----------------------------------------------------------
    
    SL.lasso.fit = ate.SL(dat$Y[1:n], dat$A[1:n], 
                          mu1.lasso$cf.pred[1:n], mu0.lasso$cf.pred[1:n], 
                          pi.lasso$cf.pred[1:n], 1-pi.lasso$cf.pred[1:n])
    
    SSL.lasso.fit = ate.SSL(dat$Y[1:n], dat$A[1:n], 
                            mu1.lasso$cf.pred, mu0.lasso$cf.pred, 
                            pi.lasso$cf.pred, 1-pi.lasso$cf.pred, 
                            imp.A.lasso$cf.pred, 
                            imp.A.lasso$cf.pred*imp.A1Y1.lasso$cf.pred,
                            (1-imp.A.lasso$cf.pred)*imp.A0Y1.lasso$cf.pred)
    
    SL.rcal.fit = ate.SL(dat$Y[1:n], dat$A[1:n],
                         mu1.rw$cf.pred[1:n], mu0.rw$cf.pred[1:n],
                         pi1.rcal$cf.pred[1:n], pi0.rcal$cf.pred[1:n])
    
    SSL.rcal.fit = ate.SSL(dat$Y[1:n], dat$A[1:n],
                           mu1.rw$cf.pred, mu0.rw$cf.pred,
                           pi1.rcal$cf.pred, pi0.rcal$cf.pred,
                           imp.A.lasso$cf.pred,
                           imp.A.lasso$cf.pred*imp.A1Y1.lasso$cf.pred,
                           (1-imp.A.lasso$cf.pred)*imp.A0Y1.lasso$cf.pred)
    
    
    nfold = nfold.ds
    foldid = c(rep_len(1:nfold, length.out = n)[rnd.n],
               rep_len(1:nfold, length.out = N-n)[rnd.N])
    
    # DS
    #-----------------------------------------------------------
    
    # Shared models
    pi.ds1 = ds.lasso(dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                      (dsid[1:n] == 1))
    pi.ds2 = ds.lasso(dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                      (dsid[1:n] == 2))
    
    # pi.bs$meas
    mu1.ds1 = ds.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                       (dsid[1:n] == 1),
                       subset = dat$A[1:n]==1)
    mu1.ds2 = ds.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                       (dsid[1:n] == 2),
                       subset = dat$A[1:n]==1)
    # mu1.bs$meas
    mu0.ds1 = ds.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                       (dsid[1:n] == 1),
                       subset = dat$A[1:n]==0)
    mu0.ds2 = ds.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                       (dsid[1:n] == 2),
                       subset = dat$A[1:n]==0)
    
    # mu0.bs$meas
    pi1.rcal.ds1 = cf.rcal(dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                           subset = (dsid[1:n]==1),
                           weight = mu1.ds2$cf.pred *(1-mu1.ds2$cf.pred),
                           dcf.weight = mu1.ds2$ds.pred*(1-mu1.ds2$ds.pred))
    pi1.rcal.ds2 = cf.rcal(dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                           subset = (dsid[1:n]==2),
                           weight = mu1.ds1$cf.pred *(1-mu1.ds1$cf.pred),
                           dcf.weight = mu1.ds1$ds.pred*(1-mu1.ds1$ds.pred))
    pi0.rcal.ds1 = cf.rcal(1-dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                           subset = (dsid[1:n]==1),
                           weight = mu0.ds2$cf.pred *(1-mu0.ds2$cf.pred),
                           dcf.weight = mu0.ds2$ds.pred*(1-mu0.ds2$ds.pred))
    pi0.rcal.ds2 = cf.rcal(1-dat$A[1:n], dat$X[1:n,], nfold, foldid[1:n],
                           subset = (dsid[1:n]==2),
                           weight = mu0.ds1$cf.pred *(1-mu0.ds1$cf.pred),
                           dcf.weight = mu0.ds1$ds.pred*(1-mu0.ds1$ds.pred))
    mu1.rw.ds1 = cf.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                          subset = (dat$A[1:n]==1)&(dsid[1:n]==1),
                          weight = (1/pi.ds2$cf.pred -1),
                          dcf.weight = (1/pi.ds2$ds.pred -1))
    mu1.rw.ds2 = cf.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                          subset = (dat$A[1:n]==1)&(dsid[1:n]==2),
                          weight = (1/pi.ds1$cf.pred -1),
                          dcf.weight = (1/pi.ds1$ds.pred -1))
    mu0.rw.ds1 = cf.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                          subset = (dat$A[1:n]==0) & (dsid[1:n]==1),
                          weight = (1/(1-pi.ds2$cf.pred) -1),
                          dcf.weight = (1/(1-pi.ds2$ds.pred) -1))
    mu0.rw.ds2 = cf.lasso(dat$Y[1:n], dat$X[1:n,], nfold, foldid[1:n],
                          subset = (dat$A[1:n]==0) & (dsid[1:n]==2),
                          weight = (1/(1-pi.ds1$cf.pred) -1),
                          dcf.weight = (1/(1-pi.ds1$ds.pred) -1))
    
    pi1 = apply(cbind(pi1.rcal.ds1$cf.pred, pi1.rcal.ds2$cf.pred), 1,
                mean, na.rm = TRUE)
    pi0 = apply(cbind(pi0.rcal.ds1$cf.pred, pi0.rcal.ds2$cf.pred), 1,
                mean, na.rm = TRUE)
    
    DS.rcal.fit = ate.SL(dat$Y[1:n], dat$A[1:n],
                         (mu1.rw.ds1$cf.pred+mu1.rw.ds2$cf.pred)/2, 
                         (mu0.rw.ds1$cf.pred+mu0.rw.ds2$cf.pred)/2,
                         pi1, pi0)
    
    # Unsupervised
    #-----------------------------------------------------------
    SA = as.numeric(dat$S[,1]>=quantile(dat$S[,1], prob = 1-mean(dat$A[1:n])))
    SY = as.numeric(dat$S[,2]>=quantile(dat$S[,2], prob = 1-mean(dat$Y[1:n])))
    
    # Shared models
    pi.ds1 = ds.lasso(SA, dat$X, nfold, foldid,
                      (dsid == 1))
    pi.ds2 = ds.lasso(SA, dat$X, nfold, foldid,
                      (dsid == 2))
    
    # pi.bs$meas
    mu1.ds1 = ds.lasso(SY, dat$X, nfold, foldid,
                       (dsid == 1),
                       subset = SA==1)
    mu1.ds2 = ds.lasso(SY, dat$X, nfold, foldid,
                       (dsid == 2),
                       subset = SA==1)
    # mu1.bs$meas
    mu0.ds1 = ds.lasso(SY, dat$X, nfold, foldid,
                       (dsid == 1),
                       subset = SA==0)
    mu0.ds2 = ds.lasso(SY, dat$X, nfold, foldid,
                       (dsid == 2),
                       subset = SA==0)
    
    # mu0.bs$meas
    pi1.rcal.ds1 = cf.rcal(SA, dat$X, nfold, foldid,
                           subset = (dsid==1),
                           weight = mu1.ds2$cf.pred *(1-mu1.ds2$cf.pred),
                           dcf.weight = mu1.ds2$ds.pred*(1-mu1.ds2$ds.pred))
    pi1.rcal.ds2 = cf.rcal(SA, dat$X, nfold, foldid,
                           subset = (dsid==2),
                           weight = mu1.ds1$cf.pred *(1-mu1.ds1$cf.pred),
                           dcf.weight = mu1.ds1$ds.pred*(1-mu1.ds1$ds.pred))
    pi0.rcal.ds1 = cf.rcal(1-SA, dat$X, nfold, foldid,
                           subset = (dsid==1),
                           weight = mu0.ds2$cf.pred *(1-mu0.ds2$cf.pred),
                           dcf.weight = mu0.ds2$ds.pred*(1-mu0.ds2$ds.pred))
    pi0.rcal.ds2 = cf.rcal(1-SA, dat$X, nfold, foldid,
                           subset = (dsid==2),
                           weight = mu0.ds1$cf.pred *(1-mu0.ds1$cf.pred),
                           dcf.weight = mu0.ds1$ds.pred*(1-mu0.ds1$ds.pred))
    mu1.rw.ds1 = cf.lasso(SY, dat$X, nfold, foldid,
                          subset = (SA==1)&(dsid==1),
                          weight = (1/pi.ds2$cf.pred -1),
                          dcf.weight = (1/pi.ds2$ds.pred -1))
    mu1.rw.ds2 = cf.lasso(SY, dat$X, nfold, foldid,
                          subset = (SA==1)&(dsid==2),
                          weight = (1/pi.ds1$cf.pred -1),
                          dcf.weight = (1/pi.ds1$ds.pred -1))
    mu0.rw.ds1 = cf.lasso(SY, dat$X, nfold, foldid,
                          subset = (SA==0) & (dsid==1),
                          weight = (1/(1-pi.ds2$cf.pred) -1),
                          dcf.weight = (1/(1-pi.ds2$ds.pred) -1))
    mu0.rw.ds2 = cf.lasso(SY, dat$X, nfold, foldid,
                          subset = (SA==0) & (dsid==2),
                          weight = (1/(1-pi.ds1$cf.pred) -1),
                          dcf.weight = (1/(1-pi.ds1$ds.pred) -1))
    pi1 = apply(cbind(pi1.rcal.ds1$cf.pred, pi1.rcal.ds2$cf.pred), 1,
                mean, na.rm = TRUE)
    pi0 = apply(cbind(pi0.rcal.ds1$cf.pred, pi0.rcal.ds2$cf.pred), 1,
                mean, na.rm = TRUE)

    UL.rcal.fit = ate.SL(SY, SA,
                         (mu1.rw.ds1$cf.pred+mu1.rw.ds2$cf.pred)/2, 
                         (mu0.rw.ds1$cf.pred+mu0.rw.ds2$cf.pred)/2,
                         pi1,  pi0)
    
    
    true.ate = mean(dat$mu1-dat$mu0)
    
    
    out = list(true = true.ate, 
                SSL.lasso.est = SSL.lasso.fit$est,
                SSL.lasso.se = SSL.lasso.fit$se, 
                SL.lasso.est = SL.lasso.fit$est,
                SL.lasso.se = SL.lasso.fit$se,
                SSL.rcal.est = SSL.rcal.fit$est,
                SSL.rcal.se = SSL.rcal.fit$se,
                SL.rcal.est = SL.rcal.fit$est,
                SL.rcal.se = SL.rcal.fit$se,
                DS.rcal.est = DS.rcal.fit$est,
                DS.rcal.se = DS.rcal.fit$se,
                UL.rcal.est = UL.rcal.fit$est,
                UL.rcal.se = UL.rcal.fit$se
                )
    save(out, file = paste0("tmpresult/hd_correct_res_rep",
                            irep, ".rda"))
    # return the result
    return(out)
    
  }

stopCluster(cl)

save(sim.hd, file = "simresult/hd_true_lasso.rda")


rm(list = objects())


