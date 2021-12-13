
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//function returns means of rows of data matrix
//X - data matrix p by n
//n - number of data points
//p - number of features
arma::colvec xbar_c(arma::mat& X, int n, int p) {
  // Initialize required parameters
  arma::colvec xbar(p);
  for(int i = 0; i < p; i++)
  {
    double sum = 0;
    for(int j = 0; j < n; j++)
    {
      sum = sum + X(i,j);
    }
    xbar(i) = sum/n;

  }
  //return xbar vector
  return(xbar);
}

//Soft thresholding function
//output - soft thresholded vector with each element softthresholded at value tau
arma::vec prox(arma::colvec& z, int n, double tau) {
  int i;
  double y;
  arma::colvec pro(n);
  for (i = 0; i < n; i++) {
    y = z(i);
    pro(i) = 0.0;
    //positive large reduced by tau to move towards zero
    if (y > tau){
      pro(i) = y - tau;
    }
    //negative large increased by tau
    else if (y < -tau){
      pro(i) = y + tau;
    }
  }
  return pro;
}


//function to update U centroid matrix given by U_i = \frac{1}{1 + n\tau}y_i +\frac{1}{1 + n\tau}\bar{x}
//output p by n matrix
arma::mat update_U(arma::mat& X, arma::mat& H, arma::mat& U, arma::mat& V,arma::colvec xbar, int n, int p,
                   arma::mat& Matrix_index1,arma::mat& Matrix_index2,arma::colvec sizes1, arma::colvec sizes2, int M1, int M2, int tau) {
  int i, j, k;
  double multiplier;
  double yi;
  int index1;
  //get the scalar value used for U calculation
  multiplier = 1/(1 + n*(tau));
  for (j = 0; j < n; j++) {
    for (i = 0; i < p; i++) {
      yi = X(i,j);
      //get l1 pairs
      if (sizes1(i) > 0){
        for (k = 0; k < sizes1(j); k++) {
          index1 = i + Matrix_index1(k,j);
          yi = yi +  H(i,index1) + (tau)*V(i,index1);
        }
      }
      //get l2 pairs
      if (sizes2(j) > 0){
        for (k = 0; k < sizes2(j); k++) {
          index1 = i + Matrix_index2(k,j);
          yi = yi -  H(i,index1) + (tau)*V(i,index1);
        }
      }
      //update each column of U
      U.col(i) = multiplier * yi + (1.0 - multiplier)*xbar;
    }
  }
  return U;
}

//function to update V centroid difference matrix given by V_l^{t+1} = prox_\frac{\gamma  w_l}{\tau}(u_{l1}^{t+1} - u_{l2}^{t+1} + \frac{H^{t}}{\tau})
//output p by nk matrix
arma::mat update_V(arma::mat& H, arma::mat& U, arma::mat& V,arma::colvec xbar, int nk, int p,double gamma,
                   arma::mat& w,arma::colvec index, double tau) {
  int  d;
  arma::colvec z(p);
  arma::colvec index_1 = index.col(0);
  arma::colvec index_2 = index.col(1);
  for (d = 0; d < nk; d++) {
    z = U.col(index_1(d) - 1) - U.col(index_2(d) - 1) - (1.0/(tau))*H.col(d);
    //call softthresholding function
    V.col(d) = prox(z,p,w(d,0)*(gamma)/(tau));
  }
  return V;
}


//function to update H lagrange multiplier matrix given by H_l^{t+1} = H_l^t + \tau(v_l^{t+1} - u_{l1}^{t+1} + u_{l2}^{t+1})
//output p by nk matrix
arma::mat H_update(arma::mat& X, arma::mat& H, arma::mat& index, arma::mat& U, arma::mat& V, int tau, int p, int q, int nK) {
  // Initialize required parameters
  int i;
  arma::colvec index_1 = index.col(0) ;
  //Rcpp::Rcout << index_1 << std::endl;
  arma::colvec index_2 = index.col(1) ;
  for (i = 0; i < nK; i++){
      H.col(i) = H.col(i) - tau * (U.col(index_1(i) - 1) - U.col(index_2(i) - 1) - V.col(i));
  }
  //return update_H matrix
  return(H);
}

//ADMM c++ algorithm
// [[Rcpp::export]]
Rcpp::List ADMM_c(arma::mat &X,arma::mat &H,arma::mat &U,arma::mat &V,arma::mat &index,arma::mat &Matrix_index1,
                  arma::mat &Matrix_index2,arma::vec &sizes1,arma::vec &sizes2,int M1,int M2,arma::vec &w,
                  int p, int n,int nk, double gamma,double tau,int num_iter=100){
  int i;
  arma::mat V_old(p, nk);
  arma::vec xbar = xbar_c(X,p,n);
  for(i = 0; i < num_iter; i++){
    V_old = V;
    //update U
    U = update_U(X,H,U,V,xbar,n,p,Matrix_index1,Matrix_index2,sizes1,sizes2,M1,M2,tau);
    //update V
    V = update_V(H,U,V,xbar,nk,p,gamma,w,index,tau);
    //update H
    H = H_update(X,H,index,U,V,tau,p,n,nk);
  }
  // Create named list
  return Rcpp::List::create(Rcpp::Named("U") = U,Rcpp::Named("V") = V,Rcpp::Named("H") = H);
}




