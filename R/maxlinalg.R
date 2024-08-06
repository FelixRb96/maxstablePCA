# --- Functions for max-linear algebra ---

#' Multiply two matrices with a matrix product that uses maxima instead of addition
#'
#' @description By calculating the entries with
#' \deqn{(A \diamond B)_{ij} = \max_{j=1,..., l} A_{il} B_{lj}}
#' for appropriate dimensions. Note that this operation is particularly
#' useful when working with multivariate exreme value distributions,
#' because, if the margins are standardized to standard Fréchet margins,
#' then the max-matrix product of a matrix A and a multivariate extreme value
#' distribution Z with standard Fréchet margins has the same margins up
#' to scaling.
#'
#' @name maxmatmul
#' @param A a non-negative array of dim n, k
#' @param B a non-negative array of dim k, l
#' @return A non netgative array of dim n, l.
#' The entries are given by the maximum of componentwise multiplication
#' of rows from A and columns from B. 
#' @export
#' @useDynLib maxstablePCA, .registration = TRUE, maxmatmulRC
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

