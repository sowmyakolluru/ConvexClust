
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


arma::uvec xbar_c(arma::mat& X, int p, int q) {
  // Initialize required parameters
  arma::uvec xbar(q);
  for(int i = 0; i < q; i++)
  {
    int sum = 0;
    for(int j = 0; j < p; j++)
    {
      sum = sum + X(i)(j);
    }
    xbar(i) = sum/p;

  }
  //return xbar vector
  return(xbar);
}
