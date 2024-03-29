% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxlinalg.R
\name{maxmatmulC}
\alias{maxmatmulC}
\title{Multiply two matrices with a faster (???) C function.}
\usage{
maxmatmulC(A, B)
}
\arguments{
\item{A}{a non-negative array of dim n, k}

\item{B}{a non-negative array of dim k, l}
}
\value{
A non netgative array of dim n, l. 
The entries are given by the maximum of componentwise multiplication
of rows from A and columns from B.
}
\description{
Max-stable distributions are preserved by taking 
maxima and scaling with non-neagtive numbers.
Thus when applying PCA for extremes it is 
natural to use max-matrix multiplication given
by calculating the entries with 
\eqn{A \cdot B = \bigvee_{j=1}^k A_{cdot j} B_{j \cdot}}.
}
\examples{
# Set up example matrices
A <- matrix(c(1,2,3,4,5,6), 2, 3)
B <- matrix(c(1,2,1,2,1,2), 3, 2)

# calling the function 
m1 <- maxmatmul(A, B)

# can be used for matrix-vector multiplication as well
v <- c(7,4,7)
m2 <- maxmatmul(A, v)
m3 <- maxmatmul(v,v)
}
