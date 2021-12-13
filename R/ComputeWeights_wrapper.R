#' Compute Kernel Weights
#'
#' This function computes weights w_ij using exp(−mu * ||ai − aj||_2^2). (i,j) are the pairs possible for given data matrix.
#' The weight matrix has dimensions of one column with rows of Number of pairs.
#'
#' @param X A data matrix of dimensions p x n
#' @param mu The non-negative parameter that controls the scaling of kernel weights
#'
#' @return
#' Returns a matrix of dimensions number of pairs by 1
#' @export
#'
#' @examples
#' #Test on multidimensional case
#' X = matrix(c(rep(1, 10), rep(0, 10)))
#' multipleDimensionX = t(cbind(X, X, X))
#' w <- ComputeWeights(multipleDimensionX, mu = 1)
#'
#'
ComputeWeights <- function(X, mu){

  # Call C++ ComputeWeights_c function to calculate the weights required for the algorithm
  w = ComputeWeights_c(X, mu)

  # Return the weight matrix
  return(w)
}
