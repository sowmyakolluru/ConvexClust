# Header for Rcpp and RcppArmadillo
#library(Rcpp)
#ibrary(RcppArmadillo)

# Source your C++ funcitons
#sourceCpp("src/ComputeWeights.cpp")
#sourceCpp("src/PathPreProcessing.cpp")
#sourceCpp("src/ADMM.cpp")

# Source your R wrappers
#source("R/ComputeWeights_wrapper.R")

#library(MASS)
#data(mammals)
#X <- as.matrix(mammals[,-1])
#make it p by n data matrix
#X <- t(scale(X,center=TRUE,scale=FALSE))
#number of data points in X
#n <- ncol(X)
#number of features in X
#p <- nrow(X)
#get weights matrix for the data matrix
#w = ComputeWeights_c(X, mu = 1)
#after calculating weights get call edges_c function to get index matrix, Matrix_index1 and Matrix_index2
#matrices and sizes1 and sizes2 vectors
#out1 = edges_c(w,n)
#testing xbar_c function from ADMM.cpp. This function returns means of rows of the data matrix
#xbar = xbar_c(multipleDimensionX, n, p)
#calculate primal_residual
#H = H_update(multipleDimensionX,H,index,U,V,tau,p,n,nk)
#U = update_U(multipleDimensionX, H, U, V,xbar, n, p,
          #Matrix_index1,Matrix_index2,sizes1, sizes2, M1, M2, tau)
#V = update_V(H, U, V,xbar, nk,p, gamma,
       #w,index, tau)





