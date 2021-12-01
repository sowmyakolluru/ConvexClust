
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//Function to compute weights using exp(−mu * ||ai − aj||_2^2)
// X - n x p matrix
// mu - given positive constant
// [[Rcpp::export]]
arma::mat ComputeWeights_c(const arma::mat& X, double mu) {
  // Initialize some parameters
  int n = X.n_rows;
  //calculate number of pairs possible to calculate weights matrix w_ij
  int num_pairs = n*(n-1)/2;
  //initialize weight matrix of num_pairs x 1 dimensions
  arma::mat weights(num_pairs,1, arma::fill::zeros);
  int k = 1;
  for(int i = 0; i < (n-1); i++)
  {
    for(int j = i+1; j < n; j++)
    {
      weights.row(k-1) = exp(-mu*pow(arma::norm(X.row(i) - X.row(j), "fro"),2));
      k = k + 1;
    }
  }
  //return weight matrix
  return(weights);
}
