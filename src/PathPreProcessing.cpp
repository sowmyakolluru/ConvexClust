
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//Function that performs calculations required for cluster path processing
// i integer input
// j integer input
// p integer input

// [[Rcpp::export]]
int Prepath_1c(int i, int j, int p) {

  //return the index value
  return(p*(i-1) - i*(i-1)/2 + j - i);
}

// [[Rcpp::export]]

//Function that performs calculations required for cluster path processing
// k is the vector of indexes of elements satisfying the given condition
// n is the integer passed in the function argument

arma::Mat<int> Prepath_2c(arma::uvec k, int n) {
  // Initialize required parameters
  arma::Mat<int> index_matrix(k.n_elem,2);
  for(int i = 0; i < k.n_elem; i++)
  {
    index_matrix(i,0) = ceil(0.5 *( 2*n -1 - (sqrt(pow(2*n-1,2) - 8*k(i)))));
    index_matrix(i,1) = k(i) - n*(index_matrix(i,0) - 1) + index_matrix(i,0)*(index_matrix(i,0) - 1)/2 + index_matrix(i,0);
  }
  //return index matrix
  return(index_matrix);
}

