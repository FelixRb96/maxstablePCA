# Functions for linear algebra 
# with using maxima instead of regular addition.

#' Mulitply two matrices with the max operation
#'
#' Max-stable distributions are preserved by taking 
#' maxima and scaling with non-neagtive numbers.
#' Thus when applying PCA for extremes it is 
#' natural to use max-matrix multiplication given
#' by calculating the entries with 
#' \eqn{A \cdot B = \bigvee_{j=1}^k A_{cdot j} B_{j \cdot}}.
#'
#' @param A a non-negative array of dim n, k
#' @param B a non-negative array of dim k, l
#' @return A non netgative array of dim n, l. 
#' The entries are given by the maximum of componentwise multiplication
#' of rows from A and columns from B. 
#' @export
#' @examples
#' # Set up example matrices
#' A <- matrix(c(1,2,3,4,5,6), 2, 3)
#' B <- matrix(c(1,2,1,2,1,2), 3, 2)
#' 
#' # calling the function 
#' m1 <- maxmatmul(A, B)
#'
#' # can be used for matrix-vector multiplication as well
#' v <- c(7,4,7)
#' m2 <- maxmatmul(A, v)
#' m3 <- maxmatmul(v,v)
maxmatmul <- function(A, B) {

  # check if A and B are arrays, if yes calculate their entries
  if(is.array(A) & is.array(B)) {
    dim1 <- nrow(A)
    dim2 <- ncol(B)
    result <- matrix(0, dim1, dim2)

    for(i in 1:dim1) {
      result[i,] <- apply(B, 2, function(x) max_scalarprod(A[i,], x))
    }
    return(result)
  } else {
    if(is.array(A) & !is.array(B)) {
      return(apply(A, 1, function(x) max_scalarprod(x, B)))
    } else if(!is.array(A) & is.array(B)) {
      return(apply(B, 2, function(x) max_scalarprod(x, A)))
    } else {
      return(max_scalarprod(A, B))
    }
  }
}

max_scalarprod <- function(x, y) return(max(x * y))

#' Multiply two matrices with a faster (???) C function. 
#' 
#' Max-stable distributions are preserved by taking 
#' maxima and scaling with non-neagtive numbers.
#' Thus when applying PCA for extremes it is 
#' natural to use max-matrix multiplication given
#' by calculating the entries with 
#' \eqn{A \cdot B = \bigvee_{j=1}^k A_{cdot j} B_{j \cdot}}.
#'
#' @param A a non-negative array of dim n, k
#' @param B a non-negative array of dim k, l
#' @return A non netgative array of dim n, l. 
#' The entries are given by the maximum of componentwise multiplication
#' of rows from A and columns from B. 
#' @export
#' @useDynLib maxstablePCA maxmatmulRC
#' @examples
#' # Set up example matrices
#' A <- matrix(c(1,2,3,4,5,6), 2, 3)
#' B <- matrix(c(1,2,1,2,1,2), 3, 2)
#' 
#' # calling the function 
#' m1 <- maxmatmul(A, B)
#'
#' # can be used for matrix-vector multiplication as well
#' v <- c(7,4,7)
#' m2 <- maxmatmul(A, v)
#' m3 <- maxmatmul(v,v)
maxmatmulC <- function(A, B) {

  if(is.array(A) & is.array(B)) {
    nrA <- nrow(A)
    ncA <- ncol(A)
    ncB <- ncol(B)
  } else if(is.array(A) & is.vector(B)) {
    nrA <- nrow(A)
    ncA <- ncol(A)
    ncB <- 1
  } else if(is.vector(A) & is.array(B)) {
    nrA <- 1
    ncA <- nrow(B)
    ncB <- ncol(B)
  } else {
    nrA <- 1
    ncA <- length(A)
    ncB <- 1
  }

  res <- .Call(maxmatmulRC, as.numeric(t(A)), as.numeric(t(B)), as.integer(nrA), as.integer(ncA), as.integer(ncB))
  res <- matrix(res, nrA, ncB, byrow = T)
  return(res)
}

