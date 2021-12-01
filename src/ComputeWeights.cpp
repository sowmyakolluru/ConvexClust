
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// X - n x p matrix
// K - number of clusters
// M - K x p cluster centers (always given)
// numIter - maximal number of iterations
// [[Rcpp::export]]
arma::mat ComputeWeights_c(const arma::mat& X, double mu) {
  // All input is assumed to be correct

  // Initialize some parameters
  int n = X.n_rows;
  int num_pairs = n*(n-1)/2;
  arma::mat weights(num_pairs,1, arma::fill::zeros);
  int k = 1;
  for(int i = 0; i < (n-1); i++)
  {
    for(int j = i+1; j < n; j++)
    {
      weights.row(k-1) = exp(-mu*pow(arma::norm(X.row(i) - X.row(j), "fro"),2));
      Rcpp::Rcout << -0.5*pow(arma::norm(X.row(i) - X.row(j), "fro"),2) << std::endl;
      k = k + 1;
    }
  }
  return(weights);

}
