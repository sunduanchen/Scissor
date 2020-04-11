

####################################################
#####  Augmented and Penalized Minimization    #####
#####  APML1 (L1/L2/Laplacian+L0)              #####
#####  Penalty: L0, L1, L2, Laplacian          #####
#####  Algorithm: one-step coordinate descent  #####
####################################################

APML1=function(x, y, family=c("gaussian", "binomial", "cox"), penalty=c("Lasso","Enet", "Net"), Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL, wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), nfolds=1, foldid=NULL, inzero=TRUE, isd=FALSE, iysd=FALSE, keep.beta=FALSE, ifast=TRUE, thresh=1e-7, maxit=1e+5) {

  #fcall=match.call()
  family=match.arg(family)
  penalty=match.arg(penalty)

  if (penalty=="Net" & is.null(Omega)) {
    penalty="Enet"
    cat("Enet was performed as no input of Omega")
  }
  if (penalty %in% c("Enet","Net") & alpha==1.0) {
    penalty="Lasso"
    cat("Lasso was performed as alpha=1.0")
  }

  if (alpha!=1.0) {
    if (is.null(Omega)) {
      penalty="Enet"
    } else if (!is.null(Omega)) {
      penalty="Net"
    }
  } else {
    penalty="Lasso"
  }

  wbeta=abs(wbeta)

  fit=switch(family,
             "gaussian"=LmL0(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,wbeta,sgn,isd,iysd,keep.beta,thresh,maxit),
             "binomial"=LogL0(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,wbeta,sgn,isd,keep.beta,thresh,maxit,threshP=1e-5),
             "cox"=CoxL0(x,y,Omega,alpha,lambda,nlambda,rlambda,nfolds,foldid,inzero,wbeta,sgn,isd,keep.beta,ifast,thresh,maxit))
  fit$family=family

  #fit$call=fcall
  class(fit)="APML1"
  return(fit)
}




