
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//Function to find indexes
// i
// j
// p

// [[Rcpp::export]]
int Index1_c(int i, int j, int p) {

  //return the index value
  return(p*(i-1) - i*(i-1)/2 + j - i);
}

// [[Rcpp::export]]

//Function to find indexes
// k
// n

arma::mat Index2_c(int k, int p) {
  // Initialize some parameters
  arma::mat index_matrix(1,2);
  index_matrix(0,0) = ceil(0.5 * (2 * p-1 - sqrt(pow(2*p-1,2) - 8*k)));
  index_matrix(0,1) = k - p*(index_matrix(0,0) - 1) + index_matrix(0,0)*(index_matrix(0,0) - 1)/2 + index_matrix(0,0);
  //return index matrix
  return(index_matrix);
}
