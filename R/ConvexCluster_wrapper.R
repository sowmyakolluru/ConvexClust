
#' Title
#'
#' @param X
#' @param Lambda
#' @param ix
#' @param M1
#' @param M2
#' @param s1
#' @param s2
#' @param w
#' @param gamma
#' @param nu
#' @param num_iter
#' @param tol_abs
#' @param tol_rel
#'
#' @return
#' @export
#'
#' @examples
convexADMM <- function(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,num_iter=100,tol_abs=1e-5,tol_rel=1e-4){


  out = ADMM_c(X,Lambda,ix,M1,M2,s1,s2,w,gamma,nu,num_iter=100,tol_abs=1e-5,tol_rel=1e-4)

  # Return U, V and eta at convergence
  return(list(U=out$U,V=out$V,Lambda=out$Lambda,nu=out$nu,
                   primal=out$primal[1:out$iter],dual=out$dual[1:out$iter],tol_primal=out$tol_primal[1:out$iter],
                   tol_dual=out$tol_dual[1:out$iter],iter=out$iter))
}
