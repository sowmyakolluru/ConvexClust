#' Convex Clustering Algorithm
#'
#' This function wrapper implements ADMM algorithm written in c++ to solve the convex clustering problem.
#' Given n data points in p dimensions,the key idea behind the convex clustering model
#' is that if two observations in the data points
#' belong to the same cluster, then their corresponding centroids.
#'
#'
#' @param X A data matrix of dimensions p x n
#' @param H A matrix of Lagrange Multipliers with dimensions p x nk where nk is number of pairs of weights possible
#' @param U An optional matrix of initial cluster centroids of dimensions p by n
#' @param V A matrix of centroid differences of dimensions p by nk initialized as a NULL matrix
#' @param w A vector of kernel weights of k elements
#' @param gamma Regularization parameter
#' @param tau Non-negative tuning parameter
#' @param num_iter numIter Optional integer defining maximum number of iterations for the algorithm (default is 100)
#'
#' @return
#' Returns a list of matrices U,V and H at convergence
#' @export
#'
#' @examples
#  Test on multidimensional case
#' X = matrix(c(rep(1, 10), rep(0, 10)))
#' multipleDimensionX = t(cbind(X, X, X))
#' w <- ComputeWeights(multipleDimensionX, mu = 1)
#' n <- ncol(multipleDimensionX)
#' p <- nrow(multipleDimensionX)
#' nk <- nrow(w)
#' H <- matrix(0,p,nk)
#' V <- matrix(0, nrow = p, ncol = nk)
#' gamma <- 0.5
#' tau <- 1.0
#'
#'
convexADMM <- function(X,H,U = NULL,V,w,gamma,tau,num_iter=100){

  n = ncol(X)
  p = nrow(X)
  nk = ncol(H)
  # Check whether U and H are NULL or not.
  #If NULL, initialize U based on p by n null matrix
  #check for compatibility with X dimensions.
  if (is.null(U)) {
    U <- matrix(0, nrow = p, ncol = n)
  }
  compatibility <- function(X, U)
  {
    if(nrow(U) != nrow(X) | ncol(U) != ncol(X))
    {
      stop("Wrong compatibility of centroid cluster matrix")
    }
  }
  compatibility(X,U)
  out1 = edges_c(w, n)
  index = out1$index
  Matrix_index1 = out1$Matrix_index1
  Matrix_index2 = out1$Matrix_index2
  sizes1 = out1$sizes1
  sizes2 = out1$sizes2
  M1 = nrow(Matrix_index1)
  M2 = nrow(Matrix_index2)
  # Call C++ ADMM_c function to perfrom ADMM
  out = ADMM_c(X,H,U,V,index,Matrix_index1,Matrix_index2,sizes1,sizes2,M1,M2,w,p,n,nk,gamma,tau,num_iter=100)

  # Return U, V and eta at convergence
  return(list(U=out$U,V=out$V,H=out$H))
}
