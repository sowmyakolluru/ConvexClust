\name{ConvexClust-package}
\alias{ConvexClust-package}
\alias{ConvexClust}
\docType{package}
\title{
\packageTitle{ConvexClust}
}
\description{
\packageDescription{ConvexClust}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ConvexClust}
\packageIndices{ConvexClust}
The package implements ADMM algorithm in order to perform clustering on the given n set of data points with p features. The clustering process starts with calculating weight matrix using ComputeWeights_wrapper.R and then ConvexCluster_wrapper takes these weights as one of the inputs and computes ADMM algorithm that returns centroid matrix, centroid difference matrix and Lagrange multiplier matrix.
}
\author{
\packageAuthor{ConvexClust}

Maintainer: \packageMaintainer{ConvexClust}
}
\references{
Eric C. Chi and Kenneth Lange (2015) Splitting Methods for Convex Clustering,
\emph{Journal of Computational and Graphical Statistics, Vol. 24, No. 4 (DECEMBER 2015),
pp. 994-1013},
\doi{https://www.jstor.org/stable/24737215}
}

