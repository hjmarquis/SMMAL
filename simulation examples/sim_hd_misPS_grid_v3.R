
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
source("source/sim_gen_hd.R")
source("core/cf_hd.R")


# Find the parameters for different AUCs

S.gen.A = "beta.bias.gen"
S.gen.Y = "beta.bias.gen"
Y.gen = "binary.Y"
X.gen = "X.ar"

n = 500
N = 10000
p = 500
ar.coef = c(sqrt(0.75),0.5)
coef.a = c(0,0.5, 0.25,0.125,0, rep(0,p-4))
coef.b1 = c(0.2,0.5, 0.25,0.125,0, rep(0,p-4))/2
coef.b0 = -c(0.2,0.5, 0.25,0.125,0, rep(0,p-4))/2
coef2 = c(0.125, 0.25,-0.5,0, rep(0,p-4))/2
Nrep = 500
nfold = 5
nfold.ds = 3
nseq = 200


seq.ba = seq(1,7, length.out = nseq)
seq.by = seq(1,7, length.out = nseq)

grid.file = "simresult/hd_misPS_grid.rda"

if(file.exists(grid.file))
{
  load(grid.file)
}else{
  auc.ba = auc.by = rep(0,nseq)
  
  for (i in 1:nseq)
  {
    S.par.A = list(a0 = 1, a1 = seq.ba[i], b0 = seq.ba[i], b1 = 1)
    S.par.Y = list(a0 = 1, a1 = seq.by[i], b0 = seq.by[i], b1 = 1)
    
    dat = sim.gen.hd(N, p, 
                     X.gen = X.gen, 
                     par.X = list(coef = ar.coef),
                     lp.A = "lp.wrong", 
                     par.lp.A = list(coef1=coef.a,coef2 = coef2),
                     par.lp.Y1 = list(coef=coef.b1),
                     par.lp.Y0 = list(coef=coef.b0),
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
registerDoRNG(seed = 531)
if(dir.exists("~/R-4.0.1/library"))
{
  clusterEvalQ(cl, .libPaths(c("~/R-4.0.1/library",.libPaths())))
}

# fit.fun = function(formula,...){coef(glm(formula,family = binomial,...))}
ia = as.integer(Sys.getenv('ia'))
iy = as.integer(Sys.getenv('iy'))


# for(ia in 1:ngrid)
{  
  S.par.A = list(a0 = 1, a1 = grid.ba[ia], b0 = grid.ba[ia], b1 = 1)
  # for (iy in 1:ngrid)
  {
    S.par.Y = list(a0 = 1, a1 = grid.by[iy], b0 = grid.by[iy], b1 = 1)
    
    setting = paste0("auc", round(auc.grid[ia]*100),
                     round(auc.grid[iy]*100))
    replist = c(outer(rep(1:np, 3), (1:ceiling(Nrep/np)-1)*np,'+'))[1:(Nrep*3)]
    Llist = rep_len(rep(1:3,each=np), length.out = Nrep*3)
    
    sim.np <- foreach(ipar = 1:(Nrep*3),.packages = pak.list ) %dopar%
      {
        source("core/ate_np.R", local = T)
        source("source/sim_gen_hd.R", local = T)
        source("core/cf_hd.R", local = T)
        source("core/ds_hd.R", local = T)
        
        irep = replist[ipar]
        iL = Llist[ipar]
        
        data.file =  paste0("tmpdat/hd_misPS_", setting, "_dat_rep", irep, ".rda")
        res.file = paste0("tmpresult/hd_misPS_", setting, "_res_rep",
                          irep,'L', iL,".rda")
        res.file.rep = paste0("tmpresult/hd_misPS_", setting, "_res_rep",
                          irep,'L', 1:3,".rda")
        
        if(file.exists(res.file))
        {
          return(NULL)
        }
        
        if(file.exists(data.file))
        {
          load(data.file)
        }else{
          dat = sim.gen.hd(N, p, 
                           X.gen = X.gen, 
                           par.X = list(coef = ar.coef),
                           lp.A = "lp.wrong", 
                           par.lp.A = list(coef1=coef.a,coef2 = coef2),
                           par.lp.Y1 = list(coef=coef.b1),
                           par.lp.Y0 = list(coef=coef.b0),
                           S.gen.A = S.gen.A, S.par.A = S.par.A, 
                           S.gen.Y = S.gen.Y, S.par.Y = S.par.Y)
          rnd.n = sample(n)
          rnd.N = sample(N-n)
          save(dat, rnd.n, rnd.N, file = data.file)
        }
        
        
        true.ate = mean(dat$mu1-dat$mu0)
        
        if(iL==1)
        {
          W = cbind(dat$X, dat$S)
          
          foldid = c(rep_len(1:nfold, length.out = n)[rnd.n],
                     rep_len(1:nfold, length.out = N-n)[rnd.N])
          
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
          
          out = list(true = true.ate, 
                     SSL.lasso.est = SSL.lasso.fit$est,
                     SSL.lasso.se = SSL.lasso.fit$se, 
                     SL.lasso.est = SL.lasso.fit$est,
                     SL.lasso.se = SL.lasso.fit$se,
                     SSL.rcal.est = SSL.rcal.fit$est,
                     SSL.rcal.se = SSL.rcal.fit$se,
                     SL.rcal.est = SL.rcal.fit$est,
                     SL.rcal.se = SL.rcal.fit$se
          )
          
        }else if(iL == 2)
        {
          nfold = nfold.ds
          foldid = c(rep_len(1:nfold, length.out = n)[rnd.n],
                     rep_len(1:nfold, length.out = N-n)[rnd.N])
          
          dsid = c(rep_len(1:2, length.out = n)[rnd.n],
                   rep_len(1:2, length.out = N-n)[rnd.N])

          
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
          
          out = list(true = true.ate, 
                     DS.rcal.est = DS.rcal.fit$est,
                     DS.rcal.se = DS.rcal.fit$se
          )
        }else{
          
          nfold = nfold.ds
          foldid = c(rep_len(1:nfold, length.out = n)[rnd.n],
                     rep_len(1:nfold, length.out = N-n)[rnd.N])
          
          dsid = c(rep_len(1:2, length.out = n)[rnd.n],
                   rep_len(1:2, length.out = N-n)[rnd.N])

          
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
          
          out = list(true = true.ate, 
                     UL.rcal.est = UL.rcal.fit$est,
                     UL.rcal.se = UL.rcal.fit$se
          )
        }
        
        save(out, file = res.file)
        
        if(all(file.exists(res.file.rep)))
        {
          file.remove(data.file)
        }
        
        return(NULL)
      }
  }
}

stopCluster(cl)

res.file.all = paste0("tmpresult/hd_misPS_", setting, "_res_rep",
       rep(1:Nrep,each=3),'L', 1:3,".rda")
if(all(file.exists(res.file.all)))
{
  sim.np = vector("list", Nrep)
  
  for(irep in 1:Nrep)
  {
    res.file.rep = paste0("tmpresult/hd_misPS_", setting, "_res_rep",
                          irep,'L', 1:3,".rda")
    sim.np[[irep]] = list(true = NA)
    if(file.exists(res.file.rep[1]))
    {
      load(res.file.rep[1])
      sim.np[[irep]] = out
    }else{
      sim.np[[irep]] = list(true = NA,
                            SSL.lasso.est = NA,
                            SSL.lasso.se = NA, 
                            SL.lasso.est = NA,
                            SL.lasso.se = NA,
                            SSL.rcal.est = NA,
                            SSL.rcal.se = NA,
                            SL.rcal.est = NA,
                            SL.rcal.se = NA)
      
    }
    if(file.exists(res.file.rep[2]))
    {
      load(res.file.rep[2])
      sim.np[[irep]]$true = out$true
      sim.np[[irep]] = c(sim.np[[irep]], out[-1])
    }else{
      sim.np[[irep]]$DS.rcal.est = NA
      sim.np[[irep]]$DS.rcal.se = NA
    }
    if(file.exists(res.file.rep[3]))
    {
      load(res.file.rep[3])
      sim.np[[irep]]$true = out$true
      sim.np[[irep]] = c(sim.np[[irep]], out[-1])
    }else{
      sim.np[[irep]]$UL.rcal.est = NA
      sim.np[[irep]]$UL.rcal.se = NA
    }
  }
  
  save(sim.np, file = paste0("simresult/hd_misPS_",
                             setting,"_v3.rda"))
  file.remove(res.file.all)
}


q(save = "no")



