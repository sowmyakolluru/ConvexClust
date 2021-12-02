#' Title
#'
#' @param X
#' @param mu
#'
#' @return
#' @export
#'
#' @examples
ComputeWeights <- function(X, mu){

  # Call C++ ComputeWeights_c function to calculate the weights required for the algorithm
  w = ComputeWeights_c(X, mu)

  # Return the weight matrix
  return(w)
}
