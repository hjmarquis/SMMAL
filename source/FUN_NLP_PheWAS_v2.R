#' Internal function to do inverse logit transformation
#' @noRd
g.logit = function(xx){exp(xx)/(exp(xx)+1)}

#' Internal function to do logit transformation
#' @noRd
logit = function(xx){log(xx/(1-xx))}

#' Internal function to get the probabilities that match the prevalence
#' @noRd
Prob.Scale = function(pp, prev){
    logit_pp = logit(pp);
    logit_pp[logit_pp == -Inf] = min(logit_pp[logit_pp > -Inf], na.rm=TRUE);
    logit_pp[logit_pp ==  Inf] = max(logit_pp[logit_pp < Inf], na.rm=TRUE)
    cc = uniroot(function(xx){mean(g.logit(logit_pp-xx), na.rm=TRUE)-prev},interval=c(-10,10))$root
    g.logit(logit_pp-cc)
}

#' Internal function to do MAP algorithm
#' @noRd
MAP_PheWAS_JS = function(dat=NULL, vname=NULL, yes.con=FALSE, full.output=TRUE){

    vname.log = paste(vname,"_log",sep="")
    
    prior.all = NULL
    post.all = NULL
    group.all = NULL

    family = c( rep("poisson", length(vname)),rep("gaussian", length(vname.log)) )

    name.all = c(vname,vname.log)

    for(i in seq_along(name.all)){
        tmpfm = as.formula(paste(name.all[i],"~1"))
        if(yes.con){
            ii = i%%(length(vname))
            if(ii==0){ii = length(vname)}
            vname.log.out = setdiff(vname.log, vname.log[ii])
            tmpfm2 = as.formula(paste("~", paste(c("note",vname.log.out),collapse="+")))
            tmpfm2 = FLXPmultinom(tmpfm2)
            na.ind = apply(dat,1,function(x) { any(is.na(x)) } )
            dat.tmp = dat[!na.ind,]
        }else{
            tmpfm2 = FLXPconstant()
            na.ind = is.na(dat[,name.all[i]])
            dat.tmp = dat[!na.ind,]
        }
        set.seed(1)
        n.clust = 1
        iter = 0
        while(n.clust < 2 & iter < 5){
            tmpfit = flexmix(tmpfm, k = 2,
                             model = FLXMRglmfix(fixed =~note,varFix=FALSE, family=family[i]),
                             concomitant=tmpfm2,control=list(tol = 1e-8), data=dat.tmp)
            n.clust = length(unique(tmpfit@cluster))
            iter = iter+1
        }
        
        ##if(name.all[i]=="nlp_log"){browser()}
        if(n.clust>1){
            avgcount = apply(dat.tmp[,vname.log,drop=FALSE],1,mean, na.rm=TRUE)
            tmpdiff = mean(avgcount[tmpfit@cluster==2]) - mean(avgcount[tmpfit@cluster==1])
            tmpind =  as.numeric( (tmpdiff > 0) + 1 )
            qprobtmp = qnorm(posterior(tmpfit))
            qprob = qprobtmp[,tmpind]
            qprob[is.infinite(qprob)] = -1*qprobtmp[is.infinite(qprob),3-tmpind]
            ### deal with some extreme clustering results ###
            if(sum(!is.infinite(qprob))>=2){
                qprob[qprob == Inf] = max(qprob[!is.infinite(qprob)])
                qprob[qprob == -Inf] = min(qprob[!is.infinite(qprob)])
            }else if(sum(!is.infinite(qprob))==1){
                if(qprob[!is.infinite(qprob)] >= 0){
                    qprob[qprob == Inf] = qprob[!is.infinite(qprob)]
                    qprob[qprob == -Inf] = qnorm(1/nrow(dat))
                }else{
                    qprob[qprob == Inf] = qnorm(1-1/nrow(dat))
                    qprob[qprob == -Inf] = qprob[!is.infinite(qprob)]
                }
            }else{
                qprob[qprob == Inf] = qnorm(1-1/nrow(dat))
                qprob[qprob == -Inf] = qnorm(1/nrow(dat))
            }
            #############

            qpost.tmp = data.frame("qprob" = qprob,
                               "group" = cbind(2-tmpfit@cluster,tmpfit@cluster-1)[,tmpind],
                               "prior" = tmpfit@prior[tmpind])

            qpost.tmp = as.matrix(qpost.tmp)
            qpost = matrix(NA,nrow=nrow(dat), ncol=3)
            colnames(qpost) = colnames(qpost.tmp)
            qpost[!na.ind,] = qpost.tmp

        }else{
            qpost = data.frame("qprob" = rep( qnorm(1-1/nrow(dat)), nrow(dat)),
                               "group" = rep(1,nrow(dat)),
                               "prior" = 1)
            qpost = as.matrix(qpost)
            qpost[na.ind, ] = NA
        }

        prior.all = c(prior.all,qpost[1,"prior"])
        post.all = cbind(post.all,qpost[,"qprob"])
        group.all = cbind(group.all, qpost[,"group"])
    }

    colnames(post.all) = name.all
    final.score = (apply(pnorm(post.all),1,mean,na.rm=TRUE) + pnorm(apply(post.all,1,mean,na.rm=TRUE)))/2
    final.score[is.na(final.score)] = NA
    
    avgcount = apply(dat[,vname.log,drop=FALSE],1,mean, na.rm=TRUE)
    avgpost = apply(post.all[,vname.log,drop=FALSE],1,mean, na.rm=TRUE)
    avgcount = avgcount[!is.na(avgpost)]
    avgpost = avgpost[!is.na(avgpost)]
    
    cluster.ori = kmeans(cbind(avgcount,avgpost),centers=2)
    #cluster.ori = kmeans(cbind(dat[,vname.log],post.all[,vname.log]),centers=2) ### standardize?
    class.pos = 1 + 1*(cor(cluster.ori$cluster==2, avgcount)>0)
    prev.est0 = mean(cluster.ori$cluster==class.pos)
    prev.est = (mean(final.score,na.rm=TRUE)+prev.est0)/2

    final.score = Prob.Scale(final.score, prev.est)
    cut.MAP = quantile(final.score,1-prev.est,na.rm = TRUE)
    if(full.output){
        list("scores" = final.score, "cut.MAP" = cut.MAP, "scores.all" = post.all, "variables" = vname)
    }else{
        list("scores" = final.score, "cut.MAP" = cut.MAP, "variables" = vname)
    }
    
}


#' Internal function to check parameters
#' @noRd
check_para <- function(dat.icd = NULL, dat.nlp = NULL, dat.note = NULL,
                       nm.phe = NULL,
                       p.icd = 0.001, n.icd = 100, p.nlp = 0.001, n.nlp = 100,
                       wgt.icd.nlp = c("con.5", "con1", "sd")[1],
                       yes.con = c(FALSE, TRUE)[1]){

    threshold = c(p.icd, n.icd, p.nlp, n.nlp)
    if((!is.numeric(threshold)) | (length(threshold)<4)){
        stop("Provide valid filter values to run MAP!")
    }

    if(is.null(nm.phe)){
        stop("Provide names of phenotypes (nm.phe) to run MAP!")
    }

    #if(is.null(nm.ID)){
    #   stop("Provide variable name of ID (nm.ID) for merging the ICD, NLP, and note count data!")
    # }

    #if(is.null(nm.utl)){
    #    stop("Provide variable name of health utilizaitons (nm.utl) to adjust for!")
    #}

    if(is.null(dat.icd) | is.null(dat.nlp) | is.null(dat.note)){
        stop("At least one of the ICD, NLP, or note count data are missing!")
    }
    
    
    if(any(dat.note$mat==0) ){
        stop("some note count data is zero!")
    }
    
    
    if( any(is.na(dat.icd)) | any(is.na(dat.nlp)) | any(is.na(dat.note)) ){
        stop("At least one of the ICD, NLP, or note count data have NA values!")
    }
    
    if(any(!is.logical(yes.con))){
        stop("yes.con has to be TRUE or FALSE!")
    }

    if(any(!is.element(wgt.icd.nlp,c("con.5", "con1", "sd"))) ){
        stop("wgt.icd.nlp has to be one of {con.5, con1, sd}!")
    }

}

#' @title MAP algorithm
#' @description  Main function to perform MAP algorithm to calculate predicted probabilities of positive phenotype for each patient
#' based on NLP and ICD counts adjusted for healthcare utilization.
#' @param dat.icd Data containing ICD codes counts and must have
#' nm.ID indicating the patient IDs.
#' @param dat.nlp Data containing NLP counts and must have
#' nm.ID indicating the patient IDs.
#' @param dat.note Data containing health utilization count and must have
#' nm.ID indicating the patient IDs.
#' @param nm.phe A vector of phenotypes whose results are desired.
#' @param nm.ID Variable name for the patient IDs that are shared by dat.icd,
#' dat.nlp, and dat.note.
#' @param nm.utl Variable name for the health utilization count in dat.note.
#' @param p.icd Threshold in percentage to check if ICD data is good to run MAP.
#' @param n.icd Threshold in count to check if ICD data is good to run MAP
#' (ICD has to pass both p.icd and n.icd threshold to be included ).
#' @param p.nlp Threshold in percentage to check if NLP data is good to run MAP.
#' @param n.nlp Threshold in count to check if NLP data is good to run MAP.
#' (NLP has to pass both p.nlp and n.nlp threshold to be included ).
#' @param wgt.icd.nlp Weighting scheme to create a new variable using ICD and NLP.
#' @param yes.con A logical variable indicating if concomitant is desired.
#' @param yes.nlp A logical variable indicating if patients having zero ICD, but nonzero
#' NLP are set aside.
#' @return Returns a list with following objects:
#' \item{MAP}{A data.frame with each column being the probabilities for each phenotype.}
#' \item{cut.MAP}{A vector of cutoff values for the probabilities for each phenotype.}
#' @references Integrating natural language processing for high-throughput Multimodal
#' Automated Phenotyping (MAP) with Application to PheWAS. Katherine P. Liao, Jiehuan Sun,
#' Tianrun A. Cai, Nicholas Link, Chuan Hong, Jie Huang, Jennifer Huffman, Jessica Gronsbell, Yichi Zhang,
#' Yuk-Lam Ho, Jacqueline Honerlaw, Lauren Costa, Victor Castro, Vivian Gainer, Shawn Murphy,
#' Christopher J. O’Donnell, J. Michael Gaziano, Kelly Cho, Peter Szolovits, Isaac Kohane,
#' Sheng Yu, and Tianxi Cai with the VA Million Veteran Program
#' @examples
#' ## simulate data to test the algorithm
#' ## see the vignette for more examples
#' ncase = 500 ; ncontrol = 2000
#' dd = c(rep(0,ncontrol), rep(1,ncase))
#' note = rpois(ncase+ncontrol, 20)
#' ICD = rpois(ncase+ncontrol, dd*3 + 0.1 + 0.1*log(note))
#' NLP = rpois(ncase+ncontrol, dd*3 + 0.1 + 0.1*log(note))
#' dat.icd = data.frame(ID=1:length(ICD), disease = ICD)
#' dat.nlp = data.frame(ID=1:length(ICD), disease = NLP)
#' dat.note = data.frame(ID=1:length(ICD), note = note)
#' nm.phe = 'disease' ; nm.ID = 'ID' ; nm.utl = 'note';
#' res = MAP_PheWAS_main(dat.icd = dat.icd, dat.nlp = dat.nlp, dat.note = dat.note,
#' nm.phe = nm.phe, nm.ID = nm.ID, nm.utl = nm.utl)
#' boxplot(res$MAP~dd)

MAP_PheWAS_main <- function(dat.icd = NULL, dat.nlp = NULL, dat.note = NULL,
                            nm.phe = NULL,
                            p.icd = 0.001, n.icd = 100, p.nlp = 0.001, n.nlp = 100,
                            wgt.icd.nlp = c("con.5", "con1", "sd")[1],
                            yes.con = c(FALSE, TRUE)[1], yes.nlp = c(TRUE, FALSE)[1],
                            full.output = c(TRUE, FALSE)[1]){
    #print('only keeping columns in nm.phe list')
    if(length(intersect(colnames(dat.icd$mat),nm.phe)>0)){
      dat.icd$mat = dat.icd$mat[,which(colnames(dat.icd$mat) %in% nm.phe),drop=F]
    }else{
      dat.icd$mat = dat.icd$mat[,1,drop=F]
    }
    if(length(intersect(colnames(dat.nlp$mat),nm.phe)>0)){
      dat.nlp$mat = dat.nlp$mat[,which(colnames(dat.nlp$mat) %in% nm.phe),drop=F]
    }else{
      dat.nlp$mat = dat.nlp$mat[,1,drop=F]
    }
    
    ### checking the data inputs
    check_para(dat.icd, dat.nlp, dat.note,  nm.phe,
               p.icd, n.icd, p.nlp, n.nlp, wgt.icd.nlp, yes.con)
    
    Noset = setdiff(nm.phe, union(colnames(dat.icd$mat), colnames(dat.nlp$mat)))
    
    if(length(Noset)>0){
        warning(cat("The following phenotypes do not exist in any data: \n", Noset,"\n") )
    }
    nm.phe = setdiff(nm.phe, Noset)

    #ID = sort(intersect( intersect(dat.icd[,nm.ID], dat.nlp[,nm.ID]),
    #                     dat.note[,nm.ID]))
    
    ## only consider patients who have health utilization and at least ICD or NLP data ##
    ID = intersect(dat.note$ID, union(dat.icd$ID, dat.nlp$ID)) 
    
    matid = match(ID, dat.icd$ID)
    
    ### My way of doing it
    # if(any(is.na(matid))){
    #   browser()
    # }
    # mat.null = dat.icd$mat[matid,]
    ### the old way of doing it - too slow on VINCI, but handles if there are healthcare utilization people not
    ### in the code list
    mat.null = Matrix(0,nrow=length(ID),ncol=ncol(dat.icd$mat))
    mat.null[!is.na(matid),] = dat.icd$mat[matid[!is.na(matid)],]
    
    colnames(mat.null) = colnames(dat.icd$mat)
    dat.icd = mat.null
    na.icd.ind = Matrix(0, nrow=length(ID), ncol=1)
    if(any(is.na(matid))){
        na.icd.ind[is.na(matid),] = 1
    }
    
    #dat.icd = dat.icd[match(ID, dat.icd[,nm.ID]),] # may introduce NA for patients having no ICD
    #na.icd.ind = apply(dat.icd, 1, function(x){any(is.na(x)} )
    
    matid = match(ID, dat.nlp$ID)
    mat.null = Matrix(0,nrow=length(ID),ncol=ncol(dat.nlp$mat))
    mat.null[!is.na(matid),] = dat.nlp$mat[matid[!is.na(matid)],]
    colnames(mat.null) = colnames(dat.nlp$mat)
    dat.nlp = mat.null
    na.nlp.ind = Matrix(0,nrow=length(ID),ncol=1)
    if(any(is.na(matid))){
        na.nlp.ind[is.na(matid),] = 1
    }
    
    rm(mat.null)
    # dat.nlp = dat.nlp[match(ID, dat.nlp[,nm.ID]),] # may introduce NA for patients having no NLP
    # na.nlp.ind = apply(dat.nlp, 1, function(x){any(is.na(x)} )
    
    # dat.note = dat.note[match(ID, dat.note[,nm.ID]),]
    note = log(dat.note$mat[match(ID, dat.note$ID),] + 1.0)
    
    #MAPscore = NULL
    #CUTmap = NULL
    res.all = NULL
    for(pheno in nm.phe){
        cat(pheno,"\n")
        if(is.element(pheno, colnames(dat.icd))){
            icd = dat.icd[,pheno,drop=FALSE]
        }else{icd = Matrix(0, nrow=nrow(dat.icd), ncol=1)}

        if(is.element(pheno, colnames(dat.nlp))){
            nlp = dat.nlp[,pheno,drop=FALSE]
        }else{nlp = Matrix(0, nrow=nrow(dat.nlp), ncol=1) }

        # dat = data.frame(icd=icd, nlp=nlp)
        dat = cbind(icd, nlp)
        colnames(dat) = c("icd","nlp")
        
        #if(any(is.na(dat))){
        #    stop(paste("No Missing values are allowed! Please check ", pheno, "!\n"))
        #}

        #filter.cc = apply(dat, 2, sum)
        filter.cc = colSums(dat)
        filter.pp = colMeans(dat>0, na.rm=TRUE)

        icd.ind = NULL
        nlp.ind = NULL

        ### determine whether ICD and/or NLP satisfies the filtering criteria.
        if(filter.cc['icd'] > n.icd & filter.pp['icd'] > p.icd){

            if(filter.cc['nlp'] > n.nlp & filter.pp['nlp'] > p.nlp){
                if(yes.nlp){
                    icd.ind = (dat[,'icd',drop=FALSE] > 0) | (dat[,'nlp',drop=FALSE] > 0)
                }else{
                    icd.ind = (dat[,'icd',drop=FALSE] > 0)
                }
                
                na.icd.ind.tmp = as.vector(na.icd.ind[which(icd.ind)])
                na.nlp.ind.tmp = as.vector(na.nlp.ind[which(icd.ind)])
                dat00 = as.data.frame(as.matrix(dat[which(icd.ind),]))
                dat00[na.icd.ind.tmp==1,'icd'] = NA
                dat00[na.nlp.ind.tmp==1,'nlp'] = NA
                
                note00 = as.vector(note[which(icd.ind)])

                ## choose weight to create a new variable based on ICD and NLP
                switch(wgt.icd.nlp,
                       "con.5"={tmpw = c(0.5,0.5)},
                       "con1"={tmpw = c(1,1)},
                       "sd"={tmpw = c(1/sd(dat00[,'icd'],na.rm=TRUE),1/sd(dat00[,'nlp'], na.rm=TRUE ));
                       tmpw = tmpw/sum(tmpw)})
                dat00[,'icdnlp'] = round(as.matrix(dat00[,c('icd','nlp')])%*%tmpw)
                colnm = colnames(dat00)
                cat("**** Results w/ ICD+NLP ****\n")
            }else{
                icd.ind = (dat[,'icd',drop=FALSE] > 0)
                dat00 = as.data.frame(as.matrix(dat[which(icd.ind),]))
                note00 = as.vector(note[which(icd.ind)])
                colnm = "icd"
                cat("Poor NLP Mapping for", pheno ,", Results w/ ICD only \n")
            }

        }else{

            if(filter.cc['nlp'] > n.nlp & filter.pp['nlp'] > p.nlp){
                nlp.ind = dat[,'nlp',drop=FALSE] > 0
                dat00 = as.data.frame(as.matrix(dat[which(nlp.ind),]))
                note00 = as.vector(note[which(nlp.ind)])
                colnm = "nlp"
                cat("**** Results w/ NLP only ****\n")
            }else{
                colnm = NULL
                cat("Poor ICD and NLP Mapping for",pheno, ", use ICD >=1 as Classification \n")
            }
        }
        
        MAP.p00 = Matrix(0,nrow(dat),ncol=1)
        
        ### running MAP algorithm
        if(is.null(colnm)){
           
            #MAP.p00 = 1.0*(dat[,'icd'] >= 1); cut.MAP = 0.5
            MAP.p00[dat[,'icd'] > 0 ] = 1.0; cut.MAP = 0.5;
            MAP.p00[na.icd.ind[,1]==1] = NA;
            
            if(full.output){
                res.list = list("scores" = MAP.p00, "cut.MAP" = cut.MAP,
                "scores.all" = NULL, "variables" = 'ICDge1')
            }else{
                res.list = list("scores" = MAP.p00, "cut.MAP" = cut.MAP,
                "variables" = 'ICDge1')
            }
            
        }else{
            
            dat00 = as.data.frame(cbind(dat00[,colnm], log(dat00[,colnm]+1), note=note00))
            colnames(dat00) = c(colnm, paste(colnm, "_log", sep=""), "note")
            junk = MAP_PheWAS_JS(dat = dat00, vname = colnm, yes.con = yes.con, full.output=full.output)
            # cut.MAP = junk$cut.MAP
            
            if(is.null(icd.ind)){
                MAP.p00[which(na.nlp.ind==1)] = NA
                MAP.p00[which(nlp.ind)] = junk$scores
                
                if(full.output){
                    MAP.pall = Matrix(0,nrow(dat),ncol=ncol(junk$scores.all))
                    MAP.pall[which(na.nlp.ind==1),] = NA
                    MAP.pall[which(nlp.ind),] = junk$scores.all
                    res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                    "scores.all" = MAP.pall, "variables" = junk$variables)
                }else{
                    res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                    "variables" = junk$variables)
                }
                
            }else{
                
                if(length(junk$variables) > 1){
                    
                    MAP.p00[which(icd.ind)] = junk$scores
                    if(full.output){
                        MAP.pall = Matrix(0,nrow(dat),ncol=ncol(junk$scores.all))
                        MAP.pall[which(icd.ind),] = junk$scores.all
                        res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                        "scores.all" = MAP.pall, "variables" = junk$variables)
                    }else{
                        res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                        "variables" = junk$variables)
                    }
                    
                    
                }else{
                    
                    MAP.p00[which(na.icd.ind==1)] = NA
                    MAP.p00[which(icd.ind)] = junk$scores
                    if(full.output){
                        MAP.pall = Matrix(0,nrow(dat),ncol=ncol(junk$scores.all))
                        if(any(na.icd.ind==1)){
                            MAP.pall[which(na.icd.ind==1),] = NA 
                        }
                        MAP.pall[which(icd.ind),] = junk$scores.all
                        res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                        "scores.all" = MAP.pall, "variables" = junk$variables)
                    }else{
                        res.list = list("scores" = MAP.p00, "cut.MAP" = junk$cut.MAP,
                                        "variables" = junk$variables)
                    }
                    
                }
                
            }
        }
        
        #MAPscore = cbind(MAPscore, MAP.p00)
        #CUTmap = c(CUTmap, cut.MAP)
        res.all[[which(nm.phe==pheno)]] = res.list
    }
    
    #rownames(MAPscore) = ID
    #colnames(MAPscore) = nm.phe
    names(res.all) = nm.phe
    list(ID = ID, res.all = res.all)
}

