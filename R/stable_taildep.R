# TTT
#
#


#' Calculate the estimated stable tail dependence at point x for given dataset.
#'
#' The stable tail dependence is a function describing the 
#' full dependence structure of a multivariate max-stable random
#' vector with unit Frechet margins (loc = 0, scale = 1, shape = 1). 
#' It takes a similar place to the covariance for the gaussian distribution.
#'
#' The exact formula and its theoretical properties can be found in
#' Beirlant, Jan. Statistics of Extremes. (2004). Formula (9.60) is
#' the estimator.
#' 
#' @param x, array with ncol(data) entries that is non-negative.
#' @param s, double that is non-negative. Controls the cutoff of data for estimator
#' @param data, array where the rows are representing the observations and 
#' the columns are the points sampled from. The columns are assumed to 
#' have \emph{unit Frechet} marginals.
#' @return Non negative float, representing the stable tail depndence 
#' estimator evaluated at parameter x. 
#' @seealso [maxstablePCA::transform_unitfrechet()] to make sure the marginal 
#' distributions of the columns are standardized. 
#' @export
#' @examples
#' # set up example data
#' data <- matrix(c(
#'  69.9399283, 3.960401, 69.939928,
#'  0.3952988, 4.416034,  4.416034,
#'  4.9373416,  4.937342,  1.698470,
#'  0.8265868,  8.180792,  8.180792
#' ), 4, 3)
#'
#' # call function for 3 variate random vector sample
#' stable_tail_dependence(c(1,2,3), 3, data)
stable_tail_dependence <- function(x, s, data) {
  if(is.vector(x)) {
    if(length(x) != ncol(data)) stop(paste("x should be of length", ncol(data), "but is of length", length(x)))
    return(stable_tail_dependence_univ(x, s, data))
  } else{
    std_s_data <- function(x) stable_tail_dependence_univ(x, s, data)
    return(apply(x, 1, std_s_data))
  }
}

#'
stable_tail_dependence_univ <- function(x, s, data) {
  R <- rowSums(data)
  indices_large_entries <- which(R > length(R) / s)
  fn <- function(i) max(data[i, ] * x) / R[i]
  result <- sum(sapply(indices_large_entries, fn)) / s
  return(result)
}

#' Intern function for stable tail dependence with vectorization
#' 
#' The stable tail dependence is a function describing the 
#' full dependence structure of a multivariate max-stable random
#' vector with unit Frechet margins (loc = 0, scale = 1, shape = 1). 
#' It takes a similar place to the covariance for the gaussian distribution.
#'
#' The exact formula and its theoretical properties can be found in
#' Beirlant, Jan. Statistics of Extremes. (2004). Formula (9.60) is
#' the estimator.
#' 
#' @param x, array with ncol(data) entries that is non-negative.
#' @param s, double that is non-negative. Controls the cutoff of data for estimator
#' @param data, array where the rows are representing the observations and 
#' the columns are the points sampled from. The columns are assumed to 
#' have \emph{unit Frechet} marginals.
#' @return Non negative float, representing the stable tail depndence 
#' estimator evaluated at parameters x. 
#' @seealso [maxstablePCA::transform_unitfrechet()] to make sure the marginal 
#' distributions of the columns are standardized. 
#' @useDynLib  maxstablePCA stabletaildep_vecRC
#' @export
#' @examples
#' # set up example data
#' data <- matrix(c(
#'  69.9399283, 3.960401, 69.939928,
#'  0.3952988, 4.416034,  4.416034,
#'  4.9373416,  4.937342,  1.698470,
#'  0.8265868,  8.180792,  8.180792
#' ), 4, 3)
#'
#' # call function for 3 variate random vector sample
#' stable_tail_dependence(c(1,2,3), 3, data)
stable_tail_dependenceC <- function(x, s, data) {

    # checking for dimensions of data
    if(is.vector(data)) {
      nrowdat <- 1
      ncoldat <- length(data)
    } else if(is.array(data)) {
      nrowdat <- dim(data)[1]
      ncoldat <- dim(data)[2]
    } else {
      stop("Please check data types of input")
    }

  # checking dimensions of evaluation points
  if(is.vector(x)) {
    nevals <- 1

    if(length(x) != ncol(data)) stop(paste("x should be of length", ncol(data), "but is of length", length(x)))
    return(stable_tail_dependence_univC(x, s, data))
  } else if(is.array(x)){
    nevals <- dim(x)[1]

    result <- .Call(stabletaildep_vecRC, as.numeric(t(x)), as.numeric(t(data)), as.numeric(s), as.integer(nrowdat), as.integer(ncoldat), as.integer(nevals))
    return(result)

  } else {
    stop("input x is of wring data type")
  }

}


#' Intern function for stable tail dependence before vectorization
#' 
#' @param x, array with ncol(data) entries that is non-negative.
#' @param s, double that is non-negative. Controls the cutoff of data for estimator
#' @param data, array where the rows are representing the observations and 
#' the columns are the points sampled from. The columns are assumed to 
#' have \emph{unit Frechet} marginals.
#' @return Non negative float, representing the stable tail depndence 
#' estimator evaluated at parameter x. 
#'
#' @useDynLib  maxstablePCA stabletaildepRC
stable_tail_dependence_univC <- function(x, s, data) {

  if(is.vector(data)) {
    nrowdat <- 1
    ncoldat <- length(data)
  } else if(is.array(data)) {
    nrowdat <- dim(data)[1]
    ncoldat <- dim(data)[2]
  } else {
    stop("Please check data types of input")
  }

  res <- .Call(stabletaildepRC, as.numeric(x), as.numeric(t(data)), as.numeric(s), as.integer(nrowdat), as.integer(ncoldat))
  return(res)
}
