
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


arma::colvec xbar_c(arma::mat& X, int p, int q) {
  // Initialize required parameters
  arma::colvec xbar(q);
  for(int i = 0; i < q; i++)
  {
    int sum = 0;
    for(int j = 0; j < p; j++)
    {
      sum = sum + X(i,j);
    }
    xbar(i) = sum/p;

  }
  //return xbar vector
  return(xbar);
}

double primal_residual_c(arma::mat& U, arma::mat& V, arma::mat& index,int p, int q, int nK) {
  // Initialize required parameters
  arma::mat L(q,nK);
  arma::colvec index_1 = index.col(1);
  arma::colvec index_2 = index.col(2);
  for(int i = 0; i < nK; i++){
     L.col(i) = U.col(index_1(i)) - U.col(index_2(i));
  }
  double residual_primal = sqrt(arma::accu(arma::pow(L - V, 2)));
  //return residual_primal vector
  return(residual_primal);
}

arma::colvec vector_1(arma::colvec& Y, int j, int p, arma::colvec& k, int n) {
  for(int i = 0; i < n; i++)
  {
    k(i) = p * (Y(i) - 1) - Y(i) * (Y(i) - 1)/2 + j - Y(i);
  }
  //return the index value
  return(k);
}

arma::colvec vector_2(int e, arma::colvec& Y, int p, arma::colvec& k, int n) {
  for(int i = 0; i < n; i++)
  {
    k(i) = p * (e - 1) - e*(e - 1)/2 + Y(i) - e;
  }
  //return the index value
  return(k);
}



