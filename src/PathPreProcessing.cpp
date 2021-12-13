
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

int Prepath_1c(int i, int j, int p) {

  //return the index value
  return(p*(i-1) - i*(i-1)/2 + j - i);
}



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

//Function that returns index matrices
//w is the weight matrix
//n is the number of data points that are to be clustered
// [[Rcpp::export]]
Rcpp::List edges_c(const arma::mat& w, int n){
  arma::rowvec sizes1(n);
  arma::rowvec sizes2(n);
  //returns indexes of elements where w>0
  arma::uvec Y = find(w > 0) + 1;
  //calls Prepath_2c(Y,n) which returns index matrix for supplied Y index vector and n data points
  arma::Mat<int> index = Prepath_2c(Y, n);
  //initialize matrix M1 with P.n_rows, n
  arma::mat M1(index.n_rows,n);
  //initialize matrix M2 with P.n_rows, n
  arma::mat M2(index.n_rows,n);
  for(int i = 0; i < n; i++)
  {
    //get indexes of elements where the matrix P 1st column has i values
    arma::uvec group1 = find(index.col(0) == i+1) +1 ;
    //store the number of elements in group1 at each iteration
    sizes1(i) = group1.n_elem;
    //Rcpp::Rcout << sizes1(i) << std::endl;
    if (sizes1(i) > 0) {
      for(int j = 0; j < group1.n_elem; j++)
      {
        //Rcpp::Rcout << "7" << std::endl;
        M1(j,i) = group1(j);
        //Rcpp::Rcout << "8" << std::endl;
      }
    }
    //Rcpp::Rcout << "2" << std::endl;
    //get indexes of elements where the matrix P 2nd column has i values
    arma::uvec group2 = find(index.col(1) == i+1) + 1;
    //store the number of elements in group2 at each iteration
    sizes2(i) = group2.n_elem;
    if (sizes2(i) > 0) {
      for(int j = 0; j < group2.n_elem; j++)
      {
        M2(j,i) = group2(j);
      }
    }
  }
  int g = max(sizes1);
  int h = max(sizes2);
  arma::mat Matrix_index1(g,n);
  arma::mat Matrix_index2(h,n);
  //
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < g; j++)
    {
      Matrix_index1(j,i) = M1(j,i);
    }
    for(int j = 0; j < h; j++)
    {
      Matrix_index2(j,i) = M2(j,i);
    }
  }

  // Create named list
  return Rcpp::List::create(Rcpp::Named("index") = index,Rcpp::Named("Matrix_index1") = Matrix_index1,
                            Rcpp::Named("Matrix_index2") = Matrix_index2,
                            Rcpp::Named("sizes1") = sizes1,Rcpp::Named("sizes2") = sizes2);
}

