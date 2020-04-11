

#######################
#####  Cox model  #####
#######################

CoxL0=function(x, y, Omega=NULL, alpha=1.0, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, inzero=TRUE, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), isd=FALSE, keep.beta=FALSE, ifast=TRUE, thresh=1e-6, maxit=1e+5) {
  # x=xi; y=yi
  # Omega=NULL; alpha=0.5; lambda=NULL; nlambda=10; rlambda=NULL; nfolds=10; foldid=NULL; ifast=TRUE; thresh=1e-6; maxit=1e+5; isd=FALSE; wbeta=rep(1,ncol(x)); sgn=rep(1,ncol(x))

  N0=nrow(x); p=ncol(x)
  ifast=as.integer(ifast)

  ### Adaptive
  aPen=ifelse(all(wbeta>0), TRUE, FALSE)

  ### Lambda path
  if (is.null(lambda)) {
    ilambda=1
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p, 0.0001, 0.01)
    }
    lambda=(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
  } else {
    ilambda=0
    nlambda=length(lambda)
  }


  if (is.null(Omega)) {
    penalty=ifelse(alpha==1, "Lasso", "Enet")
    adaptive=ifelse(any(wbeta!=1), TRUE, FALSE)
  } else {
    penalty=ifelse(alpha==1, "Lasso", "Net")
    adaptive=c(ifelse(any(wbeta!=1), TRUE, FALSE),ifelse(any(sgn!=1), TRUE, FALSE))

    Omega=rbind(0,cbind(0,Omega)); sgn1=c(1,sgn) # intercept
    ## Check Omega: positive off-diagonal, zero diagonal
    if (any(diag(Omega)!=0)) {
      diag(Omega)=0
      # cat("Diagonal of Omega was set to all zeros\n")
    }
    if (any(Omega<0)) {
      Omega=abs(Omega)
      # cat("Off-diagonal of Omega was foced to non-negative values\n")
    }

    ### Correlation/Adjacency matrix
    if (inherits(Omega, "dgCMatrix")) {
      W=OmegaSC(Omega, sgn1); W$loc=W$loc+1
    } else {
      W=OmegaC(Omega, sgn1); W$loc=W$loc+1
    }
    rm(Omega)
  }


  #####  Run  #####
  prep0=PrepCox(x, y)
  out=switch(penalty,
             "Net"=NetCoxC(prep0$x, prep0$tevent, alpha, lambda, nlambda, ilambda, wbeta, W$Omega, W$loc, W$nadj, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, p, N0, thresh, maxit, ifast),
             EnetCoxC(prep0$x, prep0$tevent, alpha, lambda, nlambda, ilambda, wbeta, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, p, N0, thresh, maxit, ifast)
  )
  nlambdai=out$nlambda
  if (nlambdai==0)
    return(NULL)
  lambdai=lambda[1:nlambdai]
  
  temi=out$Beta[, 1:nlambdai]; temi[is.na(temi)]=0.0
  out$Beta=Matrix(temi, sparse=TRUE)
  out$BetaSTD=Matrix(out$BetaSTD[, 1:nlambdai], sparse=TRUE)
  out$nzero=apply(out$Beta!=0, 2, sum)
  out$flag=out$flag[1:nlambdai]

  if (nfolds==1 & is.null(foldid)) {
    fit=data.frame(lambda=lambdai, nzero=out$nzero)
    if (!isd) {
      return(list(Beta=out$Beta, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    } else {
      return(list(Beta=out$BetaSTD, fit=fit, penalty=penalty, adaptive=adaptive, flag=out$flag))
    }
  } else {

    ###  Split data for cross-validation
    if (is.null(foldid)) {
      foldid=sample(rep(seq(nfolds), length=N0))
    } else {
      nfolds=max(foldid)
    }
    tb=table(foldid);N0i=numeric(nfolds)
    for (i in 1:nfolds)
      N0i[i]=sum(tb[-i])

    prepk=list()
    for (i in 1:nfolds) {
      temid=which(foldid!=i)
      prepk[[i]]=PrepCox(x[temid, ], y[temid, ])
    }
    weighti=as.vector(tapply(y[, "status"], foldid, sum))

    #####  Cross-validation estimates  #####
    outi=list(); cvPL=matrix(NA, nrow=nfolds, ncol=nlambdai)
    for (i in 1:nfolds) {
      outi[[i]]=switch(penalty,
                       "Net"=cvNetCoxC(prepk[[i]]$x, prepk[[i]]$tevent, alpha, lambdai, nlambdai, wbeta, W$Omega, W$loc, W$nadj, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, p, N0i[i], thresh, maxit, 0, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n),
                       cvEnetCoxC(prepk[[i]]$x, prepk[[i]]$tevent, alpha, lambdai, nlambdai, wbeta, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, p, N0i[i], thresh, maxit, 0, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n)
      )
      cvPL[i, 1:outi[[i]]$nlambda]=outi[[i]]$lf[1:outi[[i]]$nlambda]-outi[[i]]$ll[1:outi[[i]]$nlambda]
    }

    temi=apply(abs(cvPL)==Inf,2,sum,na.rm=TRUE)
    if (any(temi>0)) {
      nlambdai=max(min(which(temi>0))-1,1)
    }

    cvraw=cvPL/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvPL) #
    cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
    cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))

    cvraw=cvraw[1:nlambdai]
    cvm=cvm[1:nlambdai]
    cvse=cvse[1:nlambdai]

    indexi=which.max(cvm)
    indexij=which(cvm>=(cvm[indexi]-cvse[indexi]))[1]
    temi=rep("", nlambdai)
    temi[indexi]="*"# ;temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
    #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
    temCV=data.frame(lambda=lambdai[1:nlambdai], cvm=cvm, cvse=cvse, nzero=out$nzero[1:nlambdai], index=temi[1:nlambdai], stringsAsFactors=FALSE)

    if (!inzero) {
      rm(outi)
      temCV$cvm=-temCV$cvm # compatible with LM
      if (!keep.beta) {
        return(list(Beta=out$Beta[, indexi], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      } else {
        return(list(Beta=out$Beta[,1:nlambdai], fit=temCV, lambda.min=lambdai[indexi], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      }
    }
    
    #####  Cross-validation hard threshold  #####
    il0=indexi;cvm=list();cv.max=rep(NA, nlambdai)
    repeat {
      numi=out$nzero[il0]
      Betai=sapply(outi, function(x){x$Beta[, il0]})
      Betai[is.na(Betai)]=0
      BetaSTDi=sapply(outi, function(x){x$BetaSTD[, il0]})

      Betao=apply(Betai!=0, 2, sum)
      numi2=min(max(Betao), numi)

      if (numi2>0) {
        cvPL=matrix(NA, nrow=nfolds, ncol=numi2); i=1
        for (i in 1:nfolds) {
          Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
          numj=min(Betao[i], numi)
          if (numj==0) {

            cvPL[i, ]=cvTrimCoxC(c(0.0, 0.0), numj, numi2, c(0, 0), prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)

          } else {
            BetaSTDjj=BetaSTDj
            BetaSTDjj[wbeta==0]=max(abs(BetaSTDj))+1
            temo=rank(-abs(BetaSTDjj), ties.method="min")

            temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
            temo=temo[order(temo[, 1]), ]
            cvPL[i, ]=cvTrimCoxC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
          }
        }
      } else {
        cvPL=matrix(NA, nrow=nfolds, ncol=1)
        for (i in 1:nfolds)
          cvPL[i, ]=cvTrimCoxC(c(0.0, 0.0), 0, 0, c(0, 0), prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
      }

      cvraw=cvPL/weighti;nfoldi=apply(!is.na(cvraw), 2, sum);rm(cvPL) #
      cvm[[il0]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      temi=cvm[[il0]]
      if (aPen) {
        cv.max[il0]=max(temi)
      } else {
        cv.max[il0]=ifelse(length(temi)>sum(wbeta==0),max(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
      }
      # cv.max[il0]=max(cvm[[il0]])

      il1=c(il0-1, il0+1)
      for (j in 1:2) {
        if (il1[j]>=1 & il1[j]<=nlambdai) {
          if (is.na(cv.max[il1[j]])) {
            numi=out$nzero[il1[j]]
            Betai=sapply(outi, function(x){x$Beta[, il1[j]]})
            Betai[is.na(Betai)]=0
            BetaSTDi=sapply(outi, function(x){x$BetaSTD[, il1[j]]})

            Betao=apply(Betai!=0, 2, sum)
            numi2=min(max(Betao), numi)

            if (numi2>0) {
              cvPL=matrix(NA, nrow=nfolds, ncol=numi2)
              for (i in 1:nfolds ){
                Betaj=Betai[, i]; BetaSTDj=BetaSTDi[, i]
                numj=min(Betao[i], numi)
                if (numj==0) {
                  cvPL[i, ]=cvTrimCoxC(c(0.0, 0.0), numj, numi2, c(0, 0), prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
                } else {
                  BetaSTDjj=BetaSTDj
                  BetaSTDjj[wbeta==0]=max(abs(BetaSTDj))+1
                  temo=rank(-abs(BetaSTDjj), ties.method="min")

                  temo=data.frame(temo[which(temo<=numj)], which(temo<=numj))
                  temo=temo[order(temo[, 1]), ]
                  cvPL[i, ]=cvTrimCoxC(Betaj[temo[, 2]], numj, numi2, temo[, 2]-1, prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
                }
              }
            } else {
              cvPL=matrix(NA, nrow=nfolds, ncol=1)
              for(i in 1:nfolds)
                cvPL[i, ]=cvTrimCoxC(c(0.0, 0.0), 0, 0, c(0, 0), prep0$x, prep0$N, prep0$nevent, prep0$nevent1, prep0$loc1, prep0$n, prepk[[i]]$x, prepk[[i]]$N, prepk[[i]]$nevent, prepk[[i]]$nevent1, prepk[[i]]$loc1, prepk[[i]]$n, 0, 1)
            }

            cvraw=cvPL/weighti;nfoldi=apply(!is.na(cvraw), 2, sum)
            rm(cvPL)
            cvm[[il1[j]]]=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)

            temi=cvm[[il1[j]]]
            if (aPen) {
              cv.max[il1[j]]=max(temi)
            } else {
              cv.max[il1[j]]=ifelse(length(temi)>sum(wbeta==0),max(temi[-c(1:sum(wbeta==0))]),temi[sum(wbeta==0)])
            }
            # cv.max[il1[j]]=max(cvm[[il1[j]]])
          }
        } else {
          break
        }
      }
      # if (il1[j]==1 | il1[j]==nlambdai)
      #   break
      if (il0==which.max(cv.max)) {
        break
      } else {
        il0=which.max(cv.max)
      }
    }
    index0=which.max(cv.max)

    Beta0=out$Beta[,index0]
    BetaSTD0=out$BetaSTD[,index0]

    temi=cvm[[index0]]
    if (aPen) {
      cuti=which.max(temi)
    } else {
      cuti=ifelse(length(temi)>sum(wbeta==0),which.max(temi[-c(1:sum(wbeta==0))])+sum(wbeta==0),sum(wbeta==0))
    }
    # cuti=which.max(cvm[[index0]])

    Beta0j=out$BetaSTD[,index0]
    Beta0j[which(wbeta==0)]=max(abs(Beta0j))+1

    Beta0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0
    BetaSTD0[abs(Beta0j)<=sort(abs(Beta0j),TRUE)[cuti+1]]=0

    temCV0=data.frame(lambda=lambdai[index0],cvm=cv.max[index0],nzero=cuti)

    temCV$cvm=-temCV$cvm; temCV0$cvm=-temCV0$cvm # compatible with LM
    if (!keep.beta) {
      if (!isd) {
        return(list(Beta=out$Beta[, indexi], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      } else {
        return(list(Beta=out$BetaSTD[, indexi], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      }
    } else {
      if (!isd) {
        return(list(Beta=out$Beta[,1:nlambdai], Beta0=Beta0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      } else {
        return(list(Beta=out$BetaSTD[,1:nlambdai], Beta0=BetaSTD0, fit=temCV, fit0=temCV0, lambda.min=lambdai[indexi], lambda.opt=lambdai[index0], penalty=penalty, adaptive=adaptive, flag=out$flag[1:nlambdai]))
      }
    }

  }
}





##################################################
#####  Cox: Prepare data for log-likelihood  #####
##################################################

PrepCox=function(x, y){

  N0=nrow(x)
  oi=order(y[, "status"], decreasing=TRUE)
  x=x[oi, ];y=y[oi, ]
  oi=order(y[, "time"])
  x=x[oi, ];y=y[oi, ]

  ## remove the first censored cases
  i1=which(y[, "status"]==1);mi1=min(i1)-1
  if (mi1!=0) {
    x=x[-c(1:mi1), ];y=y[-c(1:mi1), ]
  }
  ty=y[, "time"];tevent=y[, "status"]
  N=nrow(x);n1=sum(y[, "status"])

  dty=duplicated(ty) # ties

  ### for calculation of log-likelihood
  if (any(dty)) {
    tevent0=tevent
    tevent0[which(dty)]=0

    ievent=cumsum(tevent0);loc1=which(tevent0==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=tapply(tevent==1, ievent, sum)
  } else {
    ievent=cumsum(tevent);loc1=which(tevent==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=rep(1, n)
  }

  return(list(x=x, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n))
}




