#' Convex Clustering Algorithm
#'
#' This function wrapper implements ADMM algorithm written in c++ to solve the convex clustering problem.
#' Given n data points in p dimensions,the key idea behind the convex clustering model
#' is that if two observations in the data points
#' belong to the same cluster, then their corresponding centroids.
#'
#'
#' @param X A data matrix of dimensions p x n
#' @param H An optional matrix of Lagrange Multipliers with dimensions p x k where k is number of non-zero weights
#' @param U An optional matrix of initial cluster centroids of dimensions p by n
#' @param w A vector of kernel weights of k elements
#' @param gamma Regularization parameter
#' @param tau Non-negative tuning parameter
#' @param num_iter numIter Optional integer defining maximum number of iterations for the algorithm (default is 100)
#' @param eps To monitor convergence of ADMM algorithm
#'
#' @return
#' Returns a list of matrices U,V and H at convergence
#' @export
#'
#' @examples
convexADMM <- function(X,H = NULL,U = NULL,w,gamma,tau,num_iter=100, eps = 0.001){

  # Check whether U and H are NULL or not.
  #If NULL, initialize U based on p by n null matrix
  #check for compatibility with X dimensions.
  if (is.null(U)) {
    U <- matrix(0, nrow = nrow(X), ncol = ncol(X))
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
  V <- matrix(0, nrow = nrow(X), ncol = ncol(H))
  index = out1$index
  Matrix_index1 = out1$Matrix_index1
  Matrix_index2 = out1$Matrix_index2
  sizes1 = out1$sizes1
  sizes2 = out1$sizes2
  M1 = nrow(Matrix_index1)
  M2 = nrow(Matrix_index2)
  # Call C++ ADMM_c function to perfrom ADMM
  out = ADMM_c(X,H,U,V,index,Matrix_index1,Matrix_index2,sizes1,sizes2,M1,M2,w,gamma,nu,num_iter=100, eps = 0.001)

  # Return U, V and eta at convergence
  return(list(U=out$U,V=out$V,H=out$H))
}
