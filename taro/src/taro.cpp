// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <stdlib.h>
#include <stdexcept>
#include <Rcpp.h>
#include <iostream>
#include <cstdio>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;



// [[Rcpp::export]]
int nzcount(arma::vec x) {
  vec y = nonzeros(x) ;
  return y.n_elem;
}


// [[Rcpp::export]]
arma::vec wpow(arma::vec x, double gamma0) {
  //uvec ind = find(x==0);
  vec y = pow(abs(x),-1.0*gamma0);
  y.elem(find(x==0)).zeros();
  y.elem(find(x==0)).fill(10000*y.max());
  return y;
}


// [[Rcpp::export]]
arma::mat sym_inv(arma::mat x) {
  arma::mat eigvec; arma::vec eigval;
  eig_sym(eigval, eigvec, x);
  uvec ind = find(abs(eigval)>1e-8);
  // cout << ind<< eigval << std::endl;
  return eigvec.cols(ind)* diagmat(1/eigval.elem(ind))* eigvec.cols(ind).t();
}



// [[Rcpp::export]]
double softThres(double x, double lambda) {
  return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}


// [[Rcpp::export]]
double scadThres(double z, double lambda, double k) {
  if (fabs(z) < lambda*(1 + (1/k)) ) {
    return(softThres(z,lambda/k));
  } else if (fabs(z) < 3.7*lambda){
    return(((k*2.7)/(k*2.7 -1))*softThres(z,3.7*lambda/(k*2.7)));
  } else {
    return(z);
  }
}


// [[Rcpp::export]]
arma::uvec mySdiff(arma::uvec x, arma::uvec y){
  // cout<<zeros<mat>(4,5);
  for (int j = 0; j < (int) y.n_elem; j++)
    x = x.elem(find(x != y(j)));
  return(x);
}








//  Linear constrained Penalized regression with elastic net penalty
//
//  #@param XY X'Y where X is covariance matrix and Y response vector
//  #@param XX  X'X where X is covariance matrix
//  #@param A linear constraint parameter A
//  #@param B linear constraint parameter
//  #@param Lambda1 lasso penalty parameter
//  #@param Lambda2 elastic net penalty parameter
//  #@param Mu linear constraint strength
//  #@param Nu multiplicative factor for linear constraint strength
//  #@param Beta0 initial value of regression coefficient vector
//  #@param control a list of parameters controling the fitting process
//  #@return estimated regression coefficient vector
// [[Rcpp::export]]
arma::vec bregpcdenet_Rcpp(arma::vec XY, arma::mat XX, arma::mat A, arma::vec B,
                           arma::vec Lambda1, arma::vec Lambda2,double Mu, double Nu,
                           arma::vec Beta0, int method, List control){

  double inTol = control["inTol"],diff=2*inTol;
  int inMaxIter = control["inMaxIter"],counter=0,nrA = A.n_rows;
  int p = XX.n_cols,jj;
  arma::mat AA = A.t()*A;
  arma::vec AB = A.t()*B,Beta(p),Betap(p),cp(nrA),cc = zeros<vec>(nrA),shptr(p);
  arma::vec dAA(p),dXX(p),tt(p);
  arma::vec sth;
  bool Aind = all(vectorise(A) == 0);
  Beta = Beta0;
  dXX = diagvec(XX);XX.diag().zeros();
  dAA = diagvec(AA);AA.diag().zeros();
  tt  = dXX+ Lambda2;
  shptr = tt+ Mu* dAA;
  if (method == 1){
    while ((diff > inTol) & (counter < inMaxIter)){
      cp = cc;
      Betap = Beta;
      // update Beta
      for(jj=0; jj<p; jj++){
        if(Aind){
          sth = XY(jj) - Beta.t()*XX.col(jj);// update it
        }else {
          sth = XY(jj) - Beta.t()*XX.col(jj) + cc.t()*A.col(jj) + Mu*AB(jj) - Mu*accu(AA.col(jj) % Beta);
        }
        Beta(jj) = softThres(as_scalar(sth), as_scalar(Lambda1(jj)))/as_scalar(shptr(jj));
      }
      //update cc
      cc = cp - Mu*(A*Beta-B);
      //update Mu
      Mu  = Nu*Mu;
      shptr = tt+ Mu* dAA ;
      counter = counter+1;
      diff = norm(Beta-Betap,2)/norm(Betap,2)  + norm(A*Beta-B,2);
    }
  }
  if (method == 2){
    shptr = dXX+ Mu* dAA;
    while ((diff > inTol) & (counter < inMaxIter)){
      cp = cc;
      Betap = Beta;
      // update Beta
      for(jj=0; jj<p; jj++){
        if(Aind){
          sth = XY(jj) - Beta.t()*XX.col(jj);// update it
        }else {
          sth = XY(jj) - Beta.t()*XX.col(jj) + cc.t()*A.col(jj) + Mu*AB(jj) - Mu*accu(AA.col(jj) % Beta);
        }
        Beta(jj) = scadThres(as_scalar(sth)/as_scalar(shptr(jj)), as_scalar(Lambda1(jj)),  as_scalar(shptr(jj)) );
      }
      //update cc
      cc = cp - Mu*(A*Beta-B);
      //update Mu
      Mu  = Nu*Mu;
      shptr = dXX+ Mu* dAA;
      counter = counter+1;
      diff = norm(Beta-Betap,2)/norm(Betap,2)  + norm(A*Beta-B,2);
    }
  }
  return(Beta);
}





//  Linear constrained Penalized regression with elastic net penalty when X is orthogonal
//
//  #@param XY X'Y where X is covariance matrix and Y response vector
//  #@param XX  X'X where X is covariance matrix
//  #@param A linear constraint parameter A
//  #@param B linear constraint parameter
//  #@param Lambda1 lasso penalty parameter
//  #@param Lambda2 elastic net penalty parameter
//  #@param Mu linear constraint strength
//  #@param Nu multiplicative factor for linear constraint strength
//  #@param Beta0 initial value of regression coefficient vector
//  #@param control a list of parameters controling the fitting process
//  #@return estimated regression coefficient vector
// [[Rcpp::export]]
arma::vec bregpcdenetdiag_Rcpp(arma::vec XY, arma::vec XX, arma::mat A, arma::vec B,
                               arma::vec Lambda1, arma::vec Lambda2,double Mu, double Nu,
                               arma::vec Beta0, int method, List control){

  double inTol = control["inTol"],diff=2*inTol;
  int inMaxIter = control["inMaxIter"],counter=0,nrA = A.n_rows;
  int p = XX.n_rows,jj;
  arma::mat AA = A.t()*A;
  arma::vec AB = A.t()*B,Beta(p),Betap(p),cp(nrA),cc = zeros<vec>(nrA),shptr(p);
  arma::vec dAA(p),dXX(p),tt(p);
  arma::vec sth;
  bool Aind = all(vectorise(A) == 0);
  Beta = Beta0;
  dXX = XX;
  dAA = diagvec(AA);AA.diag().zeros();
  tt  = dXX+ Lambda2;
  shptr = tt+ Mu* dAA ;
  if (method == 1){
    while ((diff > inTol) & (counter < inMaxIter)){
      cp = cc;
      Betap = Beta;
      // update Beta
      for(jj=0; jj<p; jj++){
        if(Aind){
          sth = XY(jj);// update it
        }else {
          sth = XY(jj)  + cc.t()*A.col(jj) + Mu*AB(jj) - Mu*sum(AA.col(jj) % Beta);// update it
        }
        Beta(jj) = softThres(as_scalar(sth), as_scalar(Lambda1(jj)))/as_scalar(shptr(jj));
      }

      //update cc
      cc = cp - Mu*(A*Beta-B);
      //update Mu
      Mu  = Nu*Mu;
      shptr = tt+ Mu* dAA ;
      counter = counter+1;
      diff = (norm(Beta-Betap,2))/(norm(Betap,2)) + norm(A*Beta-B,2);
    }
  }

  if (method == 2){
    shptr = dXX+ Mu* dAA;
    while ((diff > inTol) & (counter < inMaxIter)){
      cp = cc;
      Betap = Beta;
      // update Beta
      for(jj=0; jj<p; jj++){
        if(Aind){
          sth = XY(jj);// update it
        }else {
          sth = XY(jj)  + cc.t()*A.col(jj) + Mu*AB(jj) - Mu*sum(AA.col(jj) % Beta);// update it
        }
        // Beta(jj) = softThres(as_scalar(sth), as_scalar(Lambda1(jj)))/as_scalar(shptr(jj));
        Beta(jj) = scadThres(as_scalar(sth)/as_scalar(shptr(jj)), as_scalar(Lambda1(jj)),  as_scalar(shptr(jj)) );
      }
      //update cc
      cc = cp - Mu*(A*Beta-B);
      //update Mu
      Mu  = Nu*Mu;
      shptr = dXX+ Mu* dAA ;
      counter = counter+1;
      diff = (norm(Beta-Betap,2))/(norm(Betap,2)) + norm(A*Beta-B,2);
    }
  }
  return(Beta);
}




// [[Rcpp::export]]
arma::cube get_cube_prod(arma::mat &b, arma::mat &c, arma::mat &naind ){
  // cout<<zeros<mat>(4,5);
  arma::cube a(b.n_cols,c.n_cols,naind.n_cols);
  arma::mat d = c;
  for (int j = 0; j < (int) naind.n_cols; j++){
    a.slice(j) = b.t()*(d.each_col()%naind.col(j));
  }
  return(a);
}

// [[Rcpp::export]]
arma::cube get_cube_prod_inv(arma::mat &b, arma::mat &c, arma::mat &naind ){
  // cout<<zeros<mat>(4,5);
  arma::cube a(b.n_cols,c.n_cols,naind.n_cols);
  arma::mat d = c;
  for (int j = 0; j < (int) naind.n_cols; j++){
    a.slice(j) = inv(b.t()*(d.each_col()%naind.col(j)));
  }
  return(a);
}


// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// tarorrr_cpp(Ytr, X0, k, cIndex,
//             Ac, Bc,
//             Z0, matrix(0,p,q), control,
//             misind, naind2)
// Speparately missing data case
// Impute missing data with zero and iterate for initialization
// [[Rcpp::export]]
Rcpp::List taro_cpp_ur(arma::mat Y, arma::mat X0,
                       arma::mat Ac,  arma::vec  Bc,
                       int nlam, arma::vec cindex,
                       Rcpp::List initw,
                       double Dini,
                       arma::mat  Zini,
                       arma::mat  Uini,
                       arma::mat  Vini,
                       arma::mat ortV,
                       arma::vec Bv,
                       Rcpp::List control,
                       int misind,
                       int pl_method,
                       arma::mat naind,
                       arma::mat lamMat){
  int pt = X0.n_cols, q = Y.n_cols, i,j;
  int p = pt - (int) cindex.n_elem, ii, outMaxIter = control["outMaxIter"];
  int n = Y.n_rows;
  Rcpp::List out;
  arma::vec wu = initw["wu"], wv = initw["wv"];
  double  wd = initw["wd"],elp, outTol = control["outTol"] ;
  double lmif = control["lamMinFac"], elalpha  = control["elnetAlpha"];
  double lmaf = control["lamMaxFac"];
  double spu = control["spU"], spv = control["spV"];
  double mu = control["mu"], nu = control["nu"];
  Rcpp::List innerControl;
  innerControl["inTol"] = control["inTol"];
  innerControl["inMaxIter"] = control["inMaxIter"];

  // Prepare preliminary variable for the analysis
  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::mat X = X0.cols(cIndexC),Z = X0.cols(cIndex);
  arma::mat xy = X.t()*Y;
  arma::mat xtyt,xtxt;
  arma::mat C(pt,q), zk(pt-p,q);
  arma::vec vest,uk, uest,vk;
  double svk=0, suk=0, dk;
  arma::cube xx_na = get_cube_prod(X, X, naind);
  arma::cube zz_na = get_cube_prod_inv(Z, Z, naind);
  arma::vec lamv1,lambdaV1,lambdaV2,lamu1,lambdaU1,lambdaU2,Vkest,Ukest,XvYk,XuYk;


  // ------------ define  path vaariable
  arma::mat uklam = zeros(p,nlam+1),vklam = zeros(q,nlam+1),tem,tem1,Ctemp(pt,q);
  arma::mat objval = zeros(outMaxIter+1,nlam+1), BIClam = ones(4,nlam+1);
  arma::vec dklam = zeros(nlam+1),lselectSeq = zeros(nlam+1);
  arma::vec indlam = zeros(nlam+1),execTime = zeros(nlam+1);
  arma::cube zpath = zeros(pt-p,q,nlam+1),Ckpath= zeros(pt,q,nlam+1);
  arma::cube ETAkpath = zeros(n,q,nlam+1);
  arma::vec maxitc = zeros(nlam+1), convval= zeros(nlam+1);
  arma::vec relerror= zeros(nlam+1),time_prof = zeros(nlam+1);
  arma::mat XuXu; //  Av = ortV,




  // generate sequence of lambda
  // arma::mat lamMat = (xy-X.t()*Z*Zini)/(wu*wd*wv.t());
  arma::mat ablamMat = abs(lamMat);
  uvec indxx = find(ablamMat == ablamMat.max());indxx =indxx(0);
  double lmax = ablamMat.max(), SSE;
  double lmx = (lmaf*lmax);//elalpha;
  double lam,lamv2, lmn = lmax*lmif, df;
  arma::vec lamSeq =  exp(linspace<vec>(log(lmx),  log(lmn), nlam)), XvXv;

  SSE = log(accu(square(Y%naind))) ;
  BIClam = SSE*BIClam;

  // freopen("output.txt","w",stdout);
  // cout << lamSeq << lmax << std::endl;
  Uini.zeros();Vini.zeros();
  Uini(as_scalar(indxx)%p) = 1; Vini(as_scalar(indxx)/p) = 1;
  Vini = Vini*sign(lamMat(indxx)); Dini = 1.0;
  uk = Uini; dk = Dini; vk = Vini; zk = Zini;
  C.rows(cIndexC) = 0*uk*dk*vk.t();
  C.rows(cIndex) = zk;
  for(ii=0;ii < nlam; ii++){
    zpath.slice(ii) = zk;
    Ckpath.slice(ii) = C;
    ETAkpath.slice(ii) = Z*zk;
  }
  wall_clock timer;
  for(ii=0; ii < nlam; ii++){
    lam = lamSeq(ii);
    // if (df ==0 ){uk = Uini; dk = Dini; vk = Vini; zk = Zini;}
    // if (pl_method == 2){uk = Uini; dk = Dini; vk = Vini; zk = Zini;}
    C.rows(cIndexC) = uk*dk*vk.t();
    C.rows(cIndex) = zk;
    lamv1 = elalpha*lam*wd*wu;
    lamv2 = 2*(1-elalpha)*lam;
    lamu1 =elalpha*lam*wd*wv;
    // lambdaV1 = wv*(lamv1.t()*abs(uk));
    // lambdaV2 = lamv2*ones<vec>(q)*(uk.t()*uk);
    // lambdaU1 = wu*(lamu1.t()*abs(vk));
    // lambdaU2 = lamv2*ones<vec>(p)*(vk.t()*vk);


    timer.tic();
    for(i = 1; i < outMaxIter; i++){
      Ctemp = C;

      // U-Step
      tem1 = (Y-Z*C.rows(cIndex))%naind;
      tem = square(vk);
      XuYk = X.t()*tem1*vk;
      XuXu = zeros<mat>(p,p);
      for(j = 0; j < q; j++) XuXu += xx_na.slice(j)*tem(j);
      // XuXu = X.t()*(X.each_col()% (naind*square(vk)));
      lambdaU1 = wu*(lamu1.t()*abs(vk));
      lambdaU2 = lamv2*ones<vec>(p)*(vk.t()*vk);
      // cout <<  "XuYk " << accu(XuYk) <<std::endl;
      // cout <<  "XuXu " << accu(XuXu) <<std::endl;
      // cout <<  "Ac/dk " << accu(Ac/dk) <<std::endl;
      // cout <<  "Bc " << accu(Bc) <<std::endl;
      // cout <<  "lambdaU2 " << accu(lambdaU2) <<std::endl;
      // cout <<  "lambdaU1 " << accu(lambdaU1) <<std::endl;

      Ukest = bregpcdenet_Rcpp(XuYk, XuXu, Ac/dk, Bc,lambdaU1,
                               lambdaU2,mu, nu, uk*dk,pl_method,innerControl);
      // cout <<  "b" <<std::endl;
      // cout <<  "LAM = " << ii << "Uiter"<< i<<  suk <<std::endl;
      Ukest.elem( find(abs(Ukest) < 1e-7) ).zeros();
      suk = accu(abs(Ukest));
      // cout <<  "LAM = " << ii << "Uiter"<< i<<  suk <<std::endl;
      if (suk == 0) {
        lselectSeq(0) = lam;
        // dk = 0;uk.zeros();vk.zeros();
        // C.rows(cIndexC)  = uk*dk*vk.t();
        break;
      } else {
        dk = as_scalar(sqrt( Ukest.t()*Ukest));
        uk = (Ukest/dk);
      }
      // cout <<  "LAM_begin = " << ii << "Viter "<< i<< std::endl;

      // V-step
      tem = X*uk;
      XvXv = naind.t()*square(tem);
      XvYk = tem1.t()*tem;
      lambdaV1 = wv*(lamv1.t()*abs(uk));
      lambdaV2 = lamv2*ones<vec>(q)*(uk.t()*uk);
      // cout <<  "LAM_begin2 = " << ii << "Viter "<< i << std::endl;
      Vkest = bregpcdenetdiag_Rcpp(XvYk, XvXv, ortV/dk, Bv,
                                   lambdaV1, lambdaV2,mu, nu, vk*dk,
                                   pl_method,innerControl);
      Vkest.elem( find(abs(Vkest) < 1e-7) ).zeros();
      svk = accu(abs(Vkest));
      // cout <<  "LAM = " << ii << "Viter "<< i<<  svk <<std::endl;
      if (svk == 0) {
        lselectSeq(0) = lam;
        // dk = 0;uk.zeros();vk.zeros();
        // C.rows(cIndexC)  = uk*dk*vk.t();
        break;
      } else {
        dk = as_scalar(sqrt( Vkest.t()*Vkest));
        vk = (Vkest/dk);
      }
      C.rows(cIndexC)  = uk*dk*vk.t();

      // Beta step
      tem = Z.t()*((Y-X*C.rows(cIndexC))%naind);
      for(j = 0; j < q; j++) zk.col(j) = zz_na.slice(j)*tem.col(j);
      C.rows(cIndex)  = zk;

      objval(i,ii) = norm( C-Ctemp,"fro")/norm(Ctemp,"fro");
      if (objval(i,ii) < outTol) {convval(ii) = 1; break;}
    }
    elp = timer.toc();
    // cout <<  "LAM = " << ii << "Viter "<< i << "svk " << svk << std::endl;
    // cout <<  "LAM = " << ii << "Viter "<< i<< "suk  "<< suk << std::endl;
    //  add solution path information
    if (svk > 0)
      if( suk > 0){
        // if(nzcount(uk) ==1) {dk =0; uk.zeros();}
        // if(nzcount(vk) ==1) {dk =0; vk.zeros();}
        // ii = ii+1;
        uklam.col(ii) = uk;
        vklam.col(ii) = vk;
        dklam(ii) = dk;
        execTime(ii) = elp;
        zpath.slice(ii) = zk;
        Ckpath.slice(ii) = C;
        ETAkpath.slice(ii) = Z*zk + X*uk*dk*vk.t();
        lselectSeq(ii) = lam;
        maxitc(ii) = i;
        relerror(ii) = objval(i,ii);
        time_prof(ii) = elp;

        SSE = log(accu(square((Y-ETAkpath.slice(ii)) % naind) ));// accu(square(Y-X*cc));
        df = accu(abs(uk) < 1e-3) + accu(abs(vk) < 1e-3) -1 + (pt-p)*q;
        lselectSeq(ii) = lam;

        BIClam(0,ii) = SSE + (df*log((double) q*n))/(q*n); //BIC
        BIClam(1,ii) = SSE + 2*df*log((double) pt*q)/(q*n); //BICP
        BIClam(2,ii) = SSE + log(log((double) n*q))*df*log((double) pt*q)/(q*n); //GIC
        BIClam(3,ii) = SSE + 2/q/n*(df); //AIC
      }

      if( (i>1) && ((nzcount(uk) > (p*spu)) || (nzcount(vk) > (q*spv)))  ) {
        // cout <<  "LAM = " << ii << "Viter"<< i<< std::endl;
        break;}
  }
  // ii = ii-1;


  // Rcpp::List out;
  out["ukpath"] = uklam.cols(0,ii);
  out["vkpath"] = vklam.cols(0,ii);
  out["dkpath"] = arma::conv_to<arma::vec>::from(dklam.head(ii+1)) ;
  out["zpath"] = arma::conv_to<arma::cube>::from(zpath.head_slices(ii+1));
  out["etapath"] = arma::conv_to<arma::cube>::from(ETAkpath.head_slices(ii+1));
  out["objkval"] = objval.cols(0,ii);
  out["Ckpath"] = arma::conv_to<arma::cube>::from(Ckpath.head_slices(ii+1));
  out["ICkpath"] = BIClam.cols(0,ii);
  out["lamKpath"] = arma::conv_to<arma::vec>::from(lselectSeq.head(ii+1));
  out["ExecTimekpath"] = arma::conv_to<arma::vec>::from(execTime.head(ii+1));
  out["lamseq"] = lamSeq;
  out["nkpath"] = ii+1;
  out["maxit"] = maxitc;
  out["converge"] = convval;
  out["converge_error"] = relerror;
  out["time_prof"] = time_prof;
  return(out);
}


// In the initialization use ridge penalty estimate





// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// tarorrr_cpp(Ytr, X0, k, cIndex,
//             Ac, Bc,
//             Z0, matrix(0,p,q), control,
//             misind, naind2)
// Speparately missing data case
// Impute missing data with zero and iterate for initialization
// [[Rcpp::export]]
Rcpp::List tarorrr_cpp_ur(arma::mat Y, arma::mat X0, arma::vec cindex,
                          arma::mat Ac,  arma::vec  Bc,
                          arma::mat  Zini, arma::mat  Cini,
                          Rcpp::List control){
  bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["inMaxIter"], iter;
  int p = pt - (int) cindex.n_elem,j,inMaxIter = control["inMaxIter"];
  Rcpp::List out;
  double epsilon = control["inTol"],Mu = control["mu"],Nu = control["nu"],elp;
  double inTol = control["inTol"];
  // maxit = maxit/2; epsilon = epsilon;

  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::mat X = X0.cols(cIndexC),Z = X0.cols(cIndex);
  arma::mat zz = inv(Z.t()*Z), zy = Z.t()*Y, xz = X.t()*Z;
  arma::mat xx = X.t()*X, xy = X.t()*Y;
  arma::mat acac = Ac.t()*Ac, acbc = Ac.t()*Bc, xtyt,xtxt;
  arma::mat cc = zeros<mat>(Ac.n_rows,1), cp;
  arma::vec vest,ue = 1e-3*randu(p), uest,ve= 1e-3*randu(q),tem_suk;
  double svk, suk, de = 1;

  arma::mat lamMat = (xy-X.t()*Z*Zini);
  arma::mat ablamMat = abs(lamMat);
  uvec indxx = find(ablamMat == ablamMat.max());indxx =indxx(0);
  // ue(as_scalar(indxx)%p) = 1; ve(as_scalar(indxx)/p) = 1;


  arma::mat C = zeros<mat>(pt,q);
  C.rows(cIndex) = Zini; C.rows(cIndexC)  = (ue*de)*ve.t();

  wall_clock timer;
  arma::mat C_temp,Ut,Vt;
  arma::vec diffobj,dt;
  diffobj.zeros(maxit);

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    // cout << "Iter" << iter << std::endl;
    C_temp = C;

    // V-step:
    vest  = trans(ue.t()*(xy-xz*C.rows(cIndex)))/accu(square(X*ue));
    svk = norm( vest,2);
    de = svk;
    ve = vest/svk;
    // if(svk==0){
    //   de = 0;ue.zeros(p);ve.zeros(q);
    //   C.rows(cIndexC) = 0*C.rows(cIndexC);
    //   break;
    // } else {
    //   de = svk;
    //   ve = vest/svk;
    // }
    // U-step: through ridge regression setup
    Mu = control["mu"];Nu = control["nu"];
    cc = zeros<mat>(Ac.n_rows,1);
    for(j=0;j <(inMaxIter/5); j++ ){
      cp = cc;
      // xtxt = pinv(xx + Mu*acac + diagmat(1e-4*ones<vec>(p)));
      xtxt = pinv(xx + Mu*acac);
      xtyt = (xy - xz*C.rows(cIndex))*ve + Mu*acbc + Ac.t()*cc;
      uest = xtxt*xtyt;
      tem_suk = Ac*uest-Bc;
      //update cc
      cc = cp - Mu*tem_suk;
      //update Mu
      Mu  = Nu*Mu;
      if (norm(tem_suk) <  inTol) break;
    }
    // cout << "init " << j << " " << norm(tem_suk) << std::endl;
    suk = norm(uest,2);
    de = suk;
    ue = uest/de;
    C.rows(cIndexC) = (ue*de)*ve.t();
    // Z-step:
    C.rows(cIndex) = zz*(zy- xz.t()*C.rows(cIndexC));

    // convergence check
    diffobj(iter) = norm( C_temp-C,"fro")/norm(C_temp,"fro") +
      norm( cc-cp,"fro");
    if (diffobj(iter) < epsilon ) {
      converged = true;
      break;
    }
  }
  elp = timer.toc();
  out["C"] = C;
  out["u"] = ue;
  out["d"] = de;
  out["v"] = ve;
  out["diffobj"] = diffobj;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}









// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// tarorrr_cpp(Ytr, X0, k, cIndex,
//             Ac, Bc,
//             Z0, matrix(0,p,q), control,
//             misind, naind2)
// Speparately missing data case
// Impute missing data with zero and iterate for initialization
// [[Rcpp::export]]
Rcpp::List tarorrr_cpp_ur2(arma::mat Y, arma::mat X0, arma::vec cindex,
                          arma::mat Ac,  arma::vec  Bc,
                          arma::mat  Zini, arma::mat  Cini,
                          Rcpp::List control){
  bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["outMaxIter"], iter;
  int p = pt - (int) cindex.n_elem;//j,inMaxIter = control["inMaxIter"]
  Rcpp::List out;
  double epsilon = control["outTol"],elp; //Mu = control["mu"],Nu = control["nu"];
  // double inTol = control["inTol"];
  // maxit = maxit/2; epsilon = epsilon;

  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::mat X = X0.cols(cIndexC),Z = X0.cols(cIndex);
  arma::mat zz = inv(Z.t()*Z), zy = Z.t()*Y;
  arma::mat acac = Ac.t()*Ac, acbc = Ac.t()*Bc; // xtyt,xtxt;
  arma::mat X2 = X*(diagmat(ones(p)) - acac/as_scalar(Ac*Ac.t()));
  arma::mat xx = sym_inv(X2.t()*X2), xy = X2.t()*Y, xz = X2.t()*Z;
  // arma::mat cc = zeros<mat>(Ac.n_rows,1), cp;
  arma::vec ue = 1e-5*randu(p), uest,ve= 1e-5*randu(q); // vest
  double suk, de = 1; // svk


  arma::mat C = zeros<mat>(pt,q);
  C.rows(cIndex) = Zini; C.rows(cIndexC)  = (ue*de)*ve.t();

  wall_clock timer;
  arma::mat C_temp,Ut,Vt;
  arma::vec diffobj,dt;
  diffobj.zeros(maxit);

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    C_temp = C;
    // U-step:
    // Ut <- MASS::ginv(crossprod(Xt)) %*% crossprod(Xt,Yt) %*% Vt;
    uest = xx*((xy - xz*C.rows(cIndex))*ve);
    suk = norm(uest,2);
    de = suk;
    ue = uest/de;
    // cout << "Iter" << " b" << std::endl;
    // V-step:
    // temp <- crossprod(Yt,Xt %*% Ut)
    // svdd <- svd(temp); Vt <- svdd$u %*% t(svdd$v)
    svd(Ut, dt,Vt, (xy.t() - (xz*C.rows(cIndex)).t())*uest);
    // cout << "SVD: " << Ut.size() << std::endl;
    // cout << dt.size() << Vt.size() << std::endl;
    ve = Ut.col(0)*Vt.t();
    C.rows(cIndexC) = uest*ve.t();

    // Z-step:
    C.rows(cIndex) = zz*(zy- xz.t()*C.rows(cIndexC));

    // convergence check
    diffobj(iter) = norm( C_temp-C,"fro")/norm(C_temp,"fro");
    if (diffobj(iter) < epsilon ) {
      // cout << "Iter" << iter << " error "<<diffobj(iter) << std::endl;
      converged = true;
      break;
    }
  }
  elp = timer.toc();
  out["C"] = C;
  out["u"] = ue;
  out["d"] = de;
  out["v"] = ve;
  out["diffobj"] = diffobj;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}







// X;Y;Z;O;r; Zini, PhiIni; ndev[ for convergence critteria]
// tarorrr_cpp(Ytr, X0, k, cIndex,
//             Ac, Bc,
//             Z0, matrix(0,p,q), control,
//             misind, naind2)
// Speparately missing data case
// Impute missing data with zero and iterate for initialization
// [[Rcpp::export]]
Rcpp::List tarorrr(arma::mat Y, arma::mat X0, int rank,
                   arma::vec cindex, arma::mat Ac,  arma::vec  Bc,
                   arma::mat  Zini, //arma::mat  Cini,
                   Rcpp::List control){
  bool converged=false;
  int pt = X0.n_cols, q = Y.n_cols,  maxit = control["outMaxIter"], iter;
  int p = pt - (int) cindex.n_elem;//j,inMaxIter = control["inMaxIter"]
  Rcpp::List out;
  double epsilon = control["outTol"],elp; //Mu = control["mu"],Nu = control["nu"];
  // double inTol = control["inTol"];
  // maxit = maxit/2; epsilon = epsilon;

  arma::uvec cIndex =  arma::conv_to< uvec >::from(cindex-1);
  arma::uvec cIndexC = mySdiff(linspace<uvec>(0,pt-1,pt), cIndex);
  arma::mat X = X0.cols(cIndexC),Z = X0.cols(cIndex);
  arma::mat zz = inv(Z.t()*Z), zy = Z.t()*Y;
  arma::mat acac = Ac.t()*Ac, acbc = Ac.t()*Bc; // xtyt,xtxt;
  arma::mat X2 = X*(diagmat(ones(p)) - acac/as_scalar(Ac*Ac.t()));
  arma::mat xx = sym_inv(X2.t()*X2), xy = X2.t()*Y, xz = X2.t()*Z;
  // arma::mat cc = zeros<mat>(Ac.n_rows,1), cp;
  arma::mat ue = 1e-1*randu(p,rank), uest,ve= 1e-1*randu(q, rank); // vest
  arma::vec de = ones<vec>(rank);
  // double suk, de = 1; // svk


  arma::mat C = zeros<mat>(pt,q);
  C.rows(cIndex) = Zini; C.rows(cIndexC)  = (ue*diagmat(de))*ve.t();

  wall_clock timer;
  arma::mat C_temp,Ut,Vt;
  arma::vec diffobj,dt;
  diffobj.zeros(maxit);

  timer.tic();
  for(iter = 1; iter < maxit; iter++){
    C_temp = C;
    // U-step:
    // Ut <- MASS::ginv(crossprod(Xt)) %*% crossprod(Xt,Yt) %*% Vt;
    uest = xx*((xy - xz*C.rows(cIndex))*ve);
    // suk = norm(uest,2);
    // de = suk;
    // ue = uest/de;
    // cout << "Iter" << " b" << std::endl;
    // V-step:
    // temp <- crossprod(Yt,Xt %*% Ut)
    // svdd <- svd(temp); Vt <- svdd$u %*% t(svdd$v)
    svd(Ut, dt,Vt, (xy.t() - (xz*C.rows(cIndex)).t())*uest);
    // cout << "SVD: " << Ut.n_rows << "SVD: " << Ut.n_cols << std::endl;
    // cout << "SVD: " << Vt.n_rows << "SVD: " << Vt.n_cols << std::endl;
    // cout << "SVD: " << dt.n_elem << std::endl;
    // cout << dt.size() << Vt.size() << std::endl;
    ve = Ut.cols(0,rank-1)*Vt.t();
    C.rows(cIndexC) = uest*ve.t();

    // Z-step:
    C.rows(cIndex) = zz*(zy- xz.t()*C.rows(cIndexC));

    // convergence check
    diffobj(iter) = norm( C_temp-C,"fro")/norm(C_temp,"fro");
    if (diffobj(iter) < epsilon ) {
      // cout << "Iter" << iter << " error "<<diffobj(iter) << std::endl;
      converged = true;
      break;
    }
  }

  de = sqrt(sum(uest%uest)).t();
  elp = timer.toc();
  out["C"] = C;
  out["u"] = uest.each_row()/de.t();
  out["d"] = de;
  out["v"] = ve;
  out["diffobj"] = diffobj;
  out["ExecTimekpath"] = elp;
  out["maxit"] = iter;
  out["converge"] = 0;
  if(converged) out["converge"] = 1;
  return(out);
}


