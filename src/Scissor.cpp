// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

//typedef Eigen::SparseMatrix<double> SpMat;
//typedef Eigen::SparseMatrix<double>::InnerIterator InIterMat;

/*****  Center and standardize  *****/
// [[Rcpp::export]]
List scaleC(Eigen::MatrixXd X){
  int i, p=X.cols(), N=X.rows();
  Eigen::VectorXd mX(p), sdX(p);
  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
    X.col(i)/=sdX(i);
  }
  return List::create(Named("x")=X, Named("sd")=sdX, Named("m")=mX);
}

/*****  Omega  *****/
// [[Rcpp::export]]
List OmegaC(Eigen::MatrixXd & Omega, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);

  //Omega.diagonal().setZero();
  Eigen::SparseMatrix<double> OmegaS=Omega.sparseView();

  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }

  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }

  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}

/*****  Sparse Omega  *****/
// [[Rcpp::export]]
List OmegaSC(Eigen::SparseMatrix<double> & OmegaS, Eigen::VectorXi & sgn){
  int i, j, p=sgn.size();
  Eigen::VectorXi nadj=Eigen::VectorXi::Zero(p);
  Eigen::VectorXd ndegree=Eigen::VectorXd::Zero(p);

  for(i=0;i<p;++i){
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      ++nadj(i);
      ndegree(i)+=it.value();
    }
  }

  Eigen::MatrixXi loc=Eigen::MatrixXi::Zero(nadj.maxCoeff(), p);
  for(i=0;i<p;++i){
    j=0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(OmegaS, i);it;++it){
      loc(j++, i)=it.index();
      OmegaS.coeffRef(it.index(), i)=it.value()*sgn(i)*sgn(it.index())/sqrt(ndegree(i)*ndegree(it.index()));
    }
  }

  return(List::create(Named("nadj")=nadj, Named("loc")=loc, Named("Omega")=OmegaS));
}





/////////////////////////////////
/////   Linear Regression   /////
/////////////////////////////////

/*****  LM: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLmC(Eigen::MatrixXd X, Eigen::VectorXd y, double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i;
  double LiMax=0.0, LiMaxi=0.0;

  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(y.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}



/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLmC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                          Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double a0) {
  int i, j;
  Eigen::VectorXd RSS, xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),

  yF=yF.array()-a0;

  if(nn2>0){
    RSS.setZero(nn2); //nn2= # of part of data

      if(nn==0){
        RSS(0)=yF.squaredNorm();
      }else{
        for(i=0;i<nn;i++){
          j=loco(i); //   index of nonzero beta
          xbF+=XF.col(j)*beta(i); //
            RSS(i)=(yF-xbF).squaredNorm();
        }
      }

    if(nn2>nn && nn>0){
      for(i=nn;i<nn2;i++){RSS(i)=RSS(nn-1);}
    }

    if(nn2>nn && nn==0){
      for(i=nn;i<nn2;i++){RSS(i)=RSS(0);}
    }

  }else{
    RSS.setZero(1);
    RSS(0)=yF.squaredNorm();
  }

  return(RSS);
}


/*****  LM: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
             int p, int N0, double thresh, int maxit, double thresh2){

  int  i, j, it=0, il, iadd, ia=0;
  double lambda2, zi, obj0, obj1, b0, db0, objQi, objQj, rss0, rss1, RSS0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd RSS=Eigen::VectorXd::Zero(nlambda), RSQ=Eigen::VectorXd::Zero(nlambda);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }

  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          if(zi>lambda1(j)){
            // b0=(zi-lambda1(j))/(lambda2*wbeta(j)+1); // x*x/N=1
            b0=(zi-lambda1(j))/(lambda2+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            // b0=(zi+lambda1(j))/(lambda2*wbeta(j)+1);
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}


    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("rsq")=RSQ,
                      Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}


/*****  LM: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLmC(Eigen::MatrixXd X, Eigen::VectorXd y,
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
               int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double thresh2){

  int i, j, it=0, il, iadd, ia=0;
  double lambda2, zi, obj0, obj1, rss0, rss1, b0, db0, objQi, objQj, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda), BetaSTD=Eigen::MatrixXd::Zero(p,nlambda); // beta matrix for different lambdas
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXd RSSp(nlambda), RSS(nlambda), RSQ(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);
  double a0=0.0, my=0.0;
  Eigen::MatrixXd predY=Eigen::MatrixXd::Zero(NF, nlambda);

  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);


  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
//    X.col(i)/=sdX(i);
  }
  my=y.mean();
  y=y.array()-my;

  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha); // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0; objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0); // x*x/N=1
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss0-rss1)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(y.dot(X.col(i))/N);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    a0=my; xbF.setZero(NF);
    for(i=0;i<ia;i++){
      j=active(i);
      xbF+=XF.col(j)*Beta(j,il);
      a0-=mX(j)*Beta(j,il);
    }
    xbF=xbF.array()+a0;

    predY.col(il)=xbF;
    RSSp(il)=(yF-xbF).squaredNorm();
    a0S(il)=a0;

    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("predY")=predY, Named("a0S")=a0S,
                      Named("RSS")=RSS, Named("rsq")=RSQ, Named("RSSp")=RSSp, Named("nlambda")=il));
}



/*****  LM: Network (L1+La)  *****/
  // [[Rcpp::export]]
List NetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y, double alpha,
            Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
            Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
            int p, int N0, double thresh, int maxit, double thresh2){

  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N0);
  Eigen::VectorXd RSS=Eigen::VectorXd::Zero(nlambda), RSQ(nlambda);
  double xr, dbMax, lambdaMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }
  y=y.array()-y.mean();

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLmC(X, y, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLmC(X, y, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }

  RSS0=y.squaredNorm();
  obj0=RSS0/N0/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N0);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N0+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N0+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N0/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){
          flag(il)=2;
          goto exit;}
        if(it>=maxit){
          flag(il)=1;
          break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i);
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)=std::abs(y.dot(X.col(i))/N0+lambda2*zi2);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
    return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("rsq")=RSQ,
                        Named("RSS")=RSS, Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  LM: Network (L1+La) cross-validation *****/
// [[Rcpp::export]]
List cvNetLmC(Eigen::MatrixXd & X, Eigen::VectorXd & y,double alpha,
              Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
              Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
              int N, int p, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double thresh2){

  int i, j, ij, m, it=0, il, iadd, ia=0;
  double lambda2, zi, zi2, objQi=0.0, objQj, obj0, obj1, rss0, rss1, b0, db0, RSS0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd xbF(NF);
  Eigen::VectorXd RSS(nlambda), RSQ(nlambda), RSSp(nlambda);
  double xr, dbMax;
//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  double a0=0.0, my=0.0;
  Eigen::MatrixXd predY=Eigen::MatrixXd::Zero(NF, nlambda);

  Eigen::VectorXd a0S=Eigen::VectorXd::Zero(nlambda);

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N);
//    X.col(i)/=sdX(i);
  }
  my=y.mean();
  y=y.array()-my;

  RSS0=y.squaredNorm();
  obj0=RSS0/N/2.0; rss0=RSS0;

  for(i=0;i<p;++i){
    di(i)=std::abs(y.dot(X.col(i))/N);
  }

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0; dbMax=0.0; rss1=rss0; obj1=obj0;
        for(i=0;i<ia;++i){
          j=active(i);
          xr=y.dot(X.col(j));
          zi=xr/N+beta0(j);
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);} // Omega: w_kl/sqrt(d_k*d_l),L=I-Omega ; L=SLS (included sign of beta)
          }
          zi+=lambda2*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2); //  beta^T*L*beta
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2+1.0);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            y+=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              y+=db0*X.col(j);
            }
          }

          rss0+=db0*(db0*N+2.0*xr);
          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        obj0=rss0/N/2.0+objQj+objQi*lambda2/2.0;

        if(std::abs(rss1-rss0)<std::abs(thresh2*rss1)){flag(il)=0;break;}
        if(rss0 < RSS0*0.001){flag(il)=0;break;}
        if(dbMax<thresh){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh2*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){
          flag(il)=1; break;
          // goto exit;
        }
      }//while

    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i);
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)=std::abs(y.dot(X.col(i))/N+lambda2*zi2);
        if(di(i)>lambda1(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    RSS(il)=rss0; RSQ(il)=1.0-rss0/RSS0;

    a0=my; xbF.setZero(NF);
    for(i=0;i<ia;i++){
      j=active(i);
      xbF+=XF.col(j)*Beta(j,il);
      a0-=mX(j)*Beta(j,il);
    }
    xbF=xbF.array()+a0;
    predY.col(il)=xbF;

    RSSp(il)=(yF-xbF).squaredNorm();
    a0S(il)=a0;

    //if(RSQ(il) > 0.999) goto exit;
  }//for lambda

  exit:
  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("predY")=predY, Named("a0S")=a0S,
                      Named("RSS")=RSS, Named("RSSp")=RSSp, Named("rsq")=RSQ, Named("nlambda")=il));
}





///////////////////
/////   Cox   /////
///////////////////

/*****  Cox: Lambda path (max)  *****/
// [[Rcpp::export]]
double maxLambdaCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent, int N,
                     Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
                     int n, double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i, j, q;
  double denS=N, c1=0.0;
  Eigen::VectorXd lli(N);
  double LiMax=0.0, LiMaxi=0.0;

  for(i=0;i<n;i++){
    c1+=(nevent1(i)/denS);denS-=nevent(i);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){lli(j)=tevent(j)-c1;}
  }


  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(lli.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}



/*****  Derivatives of log-pl of eta (1st&2nd order),  ties  *****/
void dletaCm(Eigen::VectorXd& exb, Eigen::VectorXd& tevent, int& N,
               Eigen::VectorXi& nevent, Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1,
               int& n, Eigen::VectorXd& pl1, Eigen::VectorXd& pl2, int& ifast, int& itwo){
    int i, j, q, ipl2=0;
    double denSi, c1=0.0, c2=0.0;
    Eigen::VectorXd denS(n);

    if(ifast==0 || itwo==1)goto two;
    denSi=exb.sum();
    for(i=0;i<n;++i){
      c1+=(nevent1(i)/denSi);c2+=(nevent1(i)/pow(denSi, 2));
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
        denSi-=exb(j);
        pl1(j)=tevent(j)-exb(j)*c1;
        pl2(j)=exb(j)*(c1-exb(j)*c2);
        if(pl2(j)<=0.0)ipl2=1;
      }
    }
    if(ipl2==1){itwo=1;if(ifast==0){goto two;}}
    return;

    two:
      denSi=0.0;c1=0.0;c2=0.0;
    for(i=n-1;i>=0;--i){
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){denSi+=exb(j);}
      denS(i)=denSi;
    }
    for(i=0;i<n;++i){
      c1+=(nevent1(i)/denS(i));c2+=(nevent1(i)/pow(denS(i), 2));
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){
        pl1(j)=tevent(j)-exb(j)*c1;
        pl2(j)=exb(j)*(c1-exb(j)*c2);
    }
  }
}

/*****  Log-pl of eta,  ties  *****/
// [[Rcpp::export]]
double pletaCm(Eigen::VectorXd& xb, Eigen::VectorXd& exb, Eigen::VectorXi& nevent,
               Eigen::VectorXi& nevent1, Eigen::VectorXi& loc1, int& n, int& ifast, int& itwo){
  int i, j, q, iSS=0;
  double ll=0.0, SSi;
  Eigen::VectorXd SS(n);

  if(ifast==0 || itwo==1)goto two;
  SSi=exb.sum();
  for(i=0;i<n;++i){
    if(SSi<=0.0)iSS=1;
    for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){ll+=xb(j);}
    ll-=nevent1(i)*log(SSi);
    for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi-=exb(j);}
  }
  if(iSS==1){itwo=1;if(ifast==0){goto two;}}
  return(ll);

  two:
    ll=0.0;SSi=0.0;
    for(i=n-1;i>=0;--i){
      for(j=loc1(i)-1, q=0;q<nevent(i);j++, q++){SSi+=exb(j);}
      SS(i)=SSi;
    }
    for(i=0;i<n;++i){
      for(j=loc1(i)-1, q=0;q<nevent1(i);j++, q++){
        ll+=xb(j)-log(SS(i));
      }
    }
    return(ll);
}


/*****  Cox: Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimCoxC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                           Eigen::MatrixXd XF, int NF,
                           Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF,
                           Eigen::MatrixXd X, int N,
                           Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1, int n,
                           int ifast, int itwo){
  int i, j;
  double lli, lfi;
  Eigen::VectorXd cv, xb=Eigen::VectorXd::Zero(N), xbF=Eigen::VectorXd::Zero(NF);
  Eigen::VectorXd exb(N), exbF(NF);

  if(nn2>0){
    cv.setZero(nn2);

    if(nn==0){
      exb=(xb.array()).exp();
      lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
      exbF=(xbF.array()).exp();
      lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
      cv(0)=lfi-lli;

    }else{
      for(i=0;i<nn;i++){
        j=loco(i);
        xb+=X.col(j)*beta(i);exb=(xb.array()).exp();
        lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        xbF+=XF.col(j)*beta(i);exbF=(xbF.array()).exp();
        lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
        cv(i)=lfi-lli;
      }
    }

    if(nn2>nn && nn>0){
      for(i=nn;i<nn2;i++){cv(i)=cv(nn-1);}
    }

    if(nn2>nn && nn==0){
      for(i=nn;i<nn2;i++){cv(i)=cv(0);}
    }

  }else{
    cv.setZero(1);

    exb=(xb.array()).exp();
    lli=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
    exbF=(xbF.array()).exp();
    lfi=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
    cv(0)=lfi-lli;
  }

  return(cv);
}


/*****  Cox: Enet (L1+L2)  *****/
  // [[Rcpp::export]]
List EnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent,
              double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
              int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
              int n, int p, int N0, double thresh, int maxit, int ifast){

  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda),BetasSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double lambdaMax;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0;objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("lambda")=lambda, Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetCoxC(Eigen::MatrixXd X, Eigen::VectorXd tevent,
                double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
                int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
                int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF,
                int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){

  int i, j, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2, objQi, objQj;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda),BetasSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double mxi;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);

    mxi=XF.col(i).mean();
    XF.col(i)=XF.col(i).array()-mxi;
  }

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQi=0.0;objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
            objQi+=pow(b0, 2);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();

    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*Betas(j,il);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}



/*****  Cox: Network (L1+La)  *****/
  // [[Rcpp::export]]
List NetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha,
             Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::VectorXd wbeta,
             Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
             int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
             int n, int p, int N0, double thresh, int maxit, int ifast){

  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda), BetasSTD=Eigen::MatrixXd::Zero(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double lambdaMax;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
  }

  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaCoxC(X, tevent, N, nevent, nevent1, loc1, n, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2i*zi2;

          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();
  }//for lambda

  exit:
  if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("nlambda")=il));
}


/*****  Cox: Network (L1+La)  cross-validation  *****/
// [[Rcpp::export]]
List cvNetCoxC(Eigen::MatrixXd & X, Eigen::VectorXd tevent, double alpha,
               Eigen::VectorXd lambda, int nlambda, Eigen::VectorXd wbeta,
               Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
               int N, Eigen::VectorXi nevent, Eigen::VectorXi nevent1, Eigen::VectorXi loc1,
               int n, int p, int N0, double thresh, int maxit, int ifast, Eigen::MatrixXd XF,
               int NF, Eigen::VectorXi neventF, Eigen::VectorXi nevent1F, Eigen::VectorXi loc1F, int nF){

  int i, j, ij, m, it, il, iadd, ia=0, itwo=0;
  double lambda2, lambda2i, zi, zi2, objQi=0.0, objQj, obj0, obj1, ll0, ll1, b0, db0, PLi, PLi2;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Betas=Eigen::MatrixXd::Zero(p, nlambda), BetasSTD=Eigen::MatrixXd::Zero(p, nlambda);
  Eigen::VectorXd lambda1(p), lambda1i(p), locbeta(nlambda), locbetaF(nlambda);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  Eigen::VectorXd exb=Eigen::VectorXd::Constant(N, 1.0), xb=Eigen::VectorXd::Zero(N);
  Eigen::VectorXd exbF(NF), xbF(NF);
  Eigen::VectorXd pl1(N), pl2(N);
//  Eigen::VectorXd mX(p), sdX(p);
  Eigen::VectorXd mX(p);
  double mxi;

  for (i=0;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);

    mxi=XF.col(i).mean();
    XF.col(i)=XF.col(i).array()-mxi;
  }

  dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
  ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
  obj0=-ll0/N0;

  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta; lambda2=lambda(il)*(1.0-alpha);
    lambda1i=lambda1*N0; lambda2i=lambda2*N0;

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        objQj=0.0;
        for(i=0;i<ia;++i){
          j=active(i);
          PLi2=pl2.dot(X.col(j).cwiseAbs2());
          zi=beta0(j)*PLi2+pl1.dot(X.col(j));
          zi2=0.0;
          for(ij=0;ij<nadj(j);++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2i*zi2;

          if(zi>lambda1i(j)){
            b0=(zi-lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else if(zi<-lambda1i(j)){
            b0=(zi+lambda1i(j))/(lambda2i+PLi2);
            db0=beta0(j)-b0;
            objQi-=db0*(beta0(j)+b0-2.0*zi2);
            beta0(j)=b0;
            pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
            xb-=db0*X.col(j);
            objQj+=std::abs(b0)*lambda1(j);
          }else{
            b0=0.0;
            if(beta0(j)!=b0){
              db0=beta0(j)-b0;
              objQi-=db0*(beta0(j)+b0-2.0*zi2);
              beta0(j)=b0;
              pl1+=(pl2.cwiseProduct(X.col(j)))*db0;
              xb-=db0*X.col(j);
            }
          }
        }//for update

        ll1=ll0;obj1=obj0;
        exb=(xb.array()).exp();
        ll0=pletaCm(xb, exb, nevent, nevent1, loc1, n, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
        obj0=-ll0/N0+objQj+objQi*lambda2/2.0;

        if(std::abs(ll1-ll0)<std::abs(thresh*ll1)){flag(il)=0;break;}
        if(std::abs(obj1-obj0)<std::abs(thresh*obj1)){flag(il)=0;break;}
        if(obj0!=obj0){flag(il)=2;goto exit;}
        if(it>=maxit){flag(il)=1;goto exit;}

        dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
        if(ifast==1 && itwo==1)goto exit;
      }//while

    dletaCm(exb, tevent, N, nevent, nevent1, loc1, n, pl1, pl2, ifast, itwo);
    if(ifast==1 && itwo==1)goto exit;
    iadd=0;
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        PLi=pl1.dot(X.col(i));zi2=0.0;
        for(ij=0;ij<nadj(i);++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        PLi+=lambda2i*zi2;
        if(std::abs(PLi)>lambda1i(i)){
          active(ia)=i;iactive(i)=1;++ia;iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    locbeta(il)=ll0;
    BetasSTD.col(il)=beta0;
    Betas.col(il)=beta0.array();//sdX.array();

    xbF.setZero(NF);
    for(i=0;i<ia;i++){j=active(i);xbF+=XF.col(j)*Betas(j,il);}
    exbF=(xbF.array()).exp();
    locbetaF(il)=pletaCm(xbF, exbF, neventF, nevent1F, loc1F, nF, ifast, itwo);
  }//for lambda

  exit:
    if(ifast==1 && itwo==1 && il>0)--il;
  return(List::create(Named("Beta")=Betas, Named("BetaSTD")=BetasSTD, Named("flag")=flag,
                      Named("ll")=locbeta, Named("lf")=locbetaF, Named("nlambda")=il));
}





///////////////////////////////////
/////   Logistic Regression   /////
///////////////////////////////////

/*****  Log: Lambda path (max) inner product <xj,y> *****/
// [[Rcpp::export]]
double maxLambdaLogC(Eigen::MatrixXd X, Eigen::VectorXd Z,
                     double alpha, Eigen::VectorXd wbeta, int N0, int p){
  int i;
  double LiMax=0.0, LiMaxi=0.0;

  for (i=0; i<p; ++i) {
    if (wbeta(i) > 0.0) {
      LiMaxi=std::abs(Z.transpose()*X.col(i))/wbeta(i); // <xj,y>/N0
      if (LiMaxi > LiMax) {
        LiMax=LiMaxi;
      }
    }
  }

  LiMax=LiMax/N0/alpha;

  return(LiMax);
}


/*****  Used for CV trimming  *****/
// [[Rcpp::export]]
Eigen::VectorXd cvTrimLogC(Eigen::VectorXd beta, int nn, int nn2, Eigen::VectorXi loco,
                           Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP) {
  int i, j;
  Eigen::VectorXd Dev(nn2), xbF=Eigen::VectorXd::Zero(NF); // xb=Eigen::VectorXd::Zero(N),
  Eigen::ArrayXd p0(NF);

  for(i=0;i<nn;i++){
    j=loco(i); //   index of nonzero beta
    xbF+=XF.col(j)*beta(i);

    for (j=0; j<NF; ++j) {
      p0(j)=1.0/(1.0+exp(-xbF(j)));
      if (p0(j) < threshP) {
        p0(j)=threshP;
      } else if (p0(j) > (1.0-threshP)) {
        p0(j)=1.0-threshP;
      }
    }

    Dev(i)=(yF.array()*p0.log()+(1.0-yF.array())*(1.0-p0).log()).sum()*(-2.0);
  }


  if(nn2>nn && nn>0){
    for(i=nn;i<nn2;i++){Dev(i)=Dev(nn-1);}
  }

  return(Dev);
}


/*****  Log: Enet (L1+L2)  *****/
// [[Rcpp::export]]
List EnetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
              double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
              int p, int N0, double thresh, int maxit, double threshP){

  int  i, j, i2, it=0, il, iadd, ia=0;
  double zi, b0, db0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;


  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  // Lambda path
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLogC(X, Z, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLogC(X, Z, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }


          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update


        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

    iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(Z.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));

  }//for lambda

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  Log: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvEnetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
                double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
                int p, int N0, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP){

  int  i, j, i2, it=0, il, iadd, ia=0;
  double zi, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
  Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;

  Eigen::ArrayXd xbF(NF), pF0(NF), LLF(nlambda);
  double llF0;

  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

      iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){
        di(i)=std::abs(Z.dot(X.col(i))/N0);
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));


    // Predict Deviance
    xbF=XF*Beta.col(il);
    for (i2=0; i2<NF; ++i2) {
      pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
      if (pF0(i2) < threshP) {
        pF0(i2)=threshP;
      } else if (pF0(i2) > (1.0-threshP)) {
        pF0(i2)=1.0-threshP;
      }
    }
    LLF(il)=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  }//for lambda

  xbF=log(yF.mean()/(1.0-yF.mean()))*XF.col(0).array();
  for (i2=0; i2<NF; ++i2) {
    pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
    if (pF0(i2) < threshP) {
      pF0(i2)=threshP;
    } else if (pF0(i2) > (1.0-threshP)) {
      pF0(i2)=1.0-threshP;
    }
  }
  llF0=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("LLF")=LLF.head(il), Named("llF0")=llF0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}





/*****  Log: Network (L1+La)  *****/
// [[Rcpp::export]]
List NetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
             double alpha, Eigen::VectorXd lambda, int nlambda, int ilambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
             Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
             int p, int N0, double thresh, int maxit, double threshP){

  int  i, j, i2, m, ij, it=0, il, iadd, ia=0;
  double zi, zi2, b0, db0;
  double lambdaMax;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;


  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  // Lambda path
  if (ilambda == 1) {
    if (alpha > 0.0) {
      lambdaMax=maxLambdaLogC(X, Z, alpha, wbeta, N0, p);
    } else {
      lambdaMax=maxLambdaLogC(X, Z, 0.001, wbeta, N0, p);
    }
    lambda=lambda.array()*lambdaMax;
  }


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          zi2=0.0;
          for(ij=0; ij<nadj(j); ++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2(j)*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

    iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){

        di(i)=std::abs(Z.dot(X.col(i))/N0);

        zi2=0.0;
        for(ij=0; ij<nadj(i); ++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)+=lambda2(i)*zi2;

        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));

  }//for lambda

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}



/*****  Log: Enet (L1+L2) cross-validation  *****/
// [[Rcpp::export]]
List cvNetLogC(Eigen::MatrixXd X, Eigen::VectorXd y,
               double alpha, Eigen::VectorXd lambda, int nlambda, Eigen::ArrayXd wbeta, Eigen::ArrayXd wbetai,
               Eigen::SparseMatrix<double> & Omega, Eigen::MatrixXd loc, Eigen::VectorXi nadj,
               int p, int N0, double thresh, int maxit, Eigen::MatrixXd XF, Eigen::VectorXd yF, int NF, double threshP){

  int  i, j, i2, ij, m, it=0, il, iadd, ia=0;
  double zi, zi2, b0, db0;
  Eigen::VectorXd beta0=Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd Beta=Eigen::MatrixXd::Zero(p,nlambda),BetaSTD=Eigen::MatrixXd::Zero(p,nlambda);
  Eigen::VectorXd lambda1(p), lambda2(p);
  Eigen::VectorXi active=Eigen::VectorXi::Zero(p), iactive=Eigen::VectorXi::Zero(p);
  Eigen::VectorXi flag=Eigen::VectorXi::Zero(nlambda);
  double dbMax;

//  Eigen::VectorXd mX(p), sdX(p), di(p);
Eigen::VectorXd mX(p), di(p);

  Eigen::VectorXd Z(N0), W(N0);
  Eigen::ArrayXd xb(N0), p0(N0);
  Eigen::MatrixXd X2=Eigen::MatrixXd::Zero(N0,p);

  Eigen::VectorXd LL(nlambda);
  double wi2, xr, ll0;

  Eigen::ArrayXd xbF(NF), pF0(NF), LLF(nlambda);
  double llF0;

  //  Initial values
//  mX(0)=0.0; sdX(0)=1.0;
  mX(0)=0.0;
  X2.col(0)=X2.col(0).array()+1.0;

  for (i=1;i<p;++i) {
    mX(i)=X.col(i).mean();
    X.col(i)=X.col(i).array()-mX(i);
//    sdX(i)=sqrt(X.col(i).squaredNorm()/N0);
//    X.col(i)/=sdX(i);
    X2.col(i)=X.col(i).array()*X.col(i).array();
  }


  // Initial constant
  beta0(0)=log(y.mean()/(1.0-y.mean()));
  iactive(0)=1; active(0)=0; ia=1;
  wbetai(0)=0.0; wbeta(0)=0.0;

  xb=beta0(0)*X.col(0).array();
  for (i2=0; i2<N0; ++i2) {
    p0(i2)=1.0/(1.0+exp(-xb(i2)));
    if (p0(i2) < threshP) {
      p0(i2)=threshP;
    } else if (p0(i2) > (1.0-threshP)) {
      p0(i2)=1.0-threshP;
    }
  }
  Z=(y.array()-p0);
  for(i=0;i<p;++i){
    di(i)=std::abs(Z.dot(X.col(i))/N0);
  }

  ll0=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();


  for(il=0;il<nlambda;++il){
    lambda1=lambda(il)*alpha*wbeta*wbetai; lambda2=lambda(il)*(1.0-alpha)*wbetai; // lambda1:vector lambda*alpha, lambda2=lambda*(1-alpha)

    for(i=0;i<p;++i){
      if(iactive(i)==0){
        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia;
        }
      }
    }

    it=0;
    local:
      while(1){
        ++it;

        W=p0*(1.0-p0);
        Z=(y.array()-p0)/W.array();

        dbMax=0.0;
        for(i=0;i<ia;++i){
          j=active(i);

          wi2=W.dot(X2.col(j))/N0;
          xr=(Z.array()*W.array()).matrix().dot(X.col(j));
          zi=xr/N0+wi2*beta0(j);

          zi2=0.0;
          for(ij=0; ij<nadj(j); ++ij){
            m=loc(ij, j)-1;
            if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, j);}
          }
          zi+=lambda2(j)*zi2;

          if(zi>lambda1(j)){
            b0=(zi-lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else if(zi<-lambda1(j)){
            b0=(zi+lambda1(j))/(lambda2(j)+wi2);
            db0=b0-beta0(j); beta0(j)=b0;
          }else{
            b0=0.0; db0=0.0;
            if(beta0(j)!=b0){
              db0=b0-beta0(j); beta0(j)=b0;
            }
          }

          Z-=db0*X.col(j);
          xb+=db0*X.col(j).array();

          dbMax=std::max(dbMax, pow(db0, 2));
        }//for update

        for (i2=0; i2<N0; ++i2) {
          p0(i2)=1.0/(1.0+exp(-xb(i2)));
          if (p0(i2) < threshP) {
            p0(i2)=threshP;
          } else if (p0(i2) > (1.0-threshP)) {
            p0(i2)=1.0-threshP;
          }
        }

        if(dbMax<thresh){flag(il)=0; break;}
        if(it>=maxit){flag(il)=1; break;
        // goto exit;
        }
      }//while

      iadd=0;
    Z=(y.array()-p0);
    for(i=0;i<p;++i){
      if(iactive(i)==0){

        di(i)=std::abs(Z.dot(X.col(i))/N0);

        zi2=0.0;
        for(ij=0; ij<nadj(i); ++ij){
          m=loc(ij, i)-1;
          if(iactive(m)==1){zi2+=beta0(m)*Omega.coeffRef(m, i);}
        }
        di(i)+=lambda2(i)*zi2;

        if(di(i)>lambda1(i)){
          active(ia)=i; iactive(i)=1; ++ia; iadd=1;
        }
      }
    }
    if(iadd==1){goto local;}

    LL(il)=(y.array()*p0.log()+(1.0-y.array())*(1.0-p0).log()).mean();

    BetaSTD.col(il)=beta0;
    Beta.col(il)=beta0.array();//sdX.array();
    Beta(0,il)=Beta(0,il)-mX.dot(Beta.col(il));


    // Predict Deviance
    xbF=XF*Beta.col(il);
    for (i2=0; i2<NF; ++i2) {
      pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
      if (pF0(i2) < threshP) {
        pF0(i2)=threshP;
      } else if (pF0(i2) > (1.0-threshP)) {
        pF0(i2)=1.0-threshP;
      }
    }
    LLF(il)=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  }//for lambda

  xbF=log(yF.mean()/(1.0-yF.mean()))*XF.col(0).array();
  for (i2=0; i2<NF; ++i2) {
    pF0(i2)=1.0/(1.0+exp(-xbF(i2)));
    if (pF0(i2) < threshP) {
      pF0(i2)=threshP;
    } else if (pF0(i2) > (1.0-threshP)) {
      pF0(i2)=1.0-threshP;
    }
  }
  llF0=(yF.array()*pF0.log()+(1.0-yF.array())*(1.0-pF0).log()).mean();

  return(List::create(Named("Beta")=Beta, Named("BetaSTD")=BetaSTD, Named("flag")=flag, Named("it")=it,
                      Named("LL")=LL.head(il), Named("ll0")=ll0,
                      Named("LLF")=LLF.head(il), Named("llF0")=llF0,
                      Named("lambda")=lambda, Named("nlambda")=il));
}




