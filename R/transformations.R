# Tools for transforming data to the necessary margins
# and given an analyzed transformed dataset the possiblity
# to transform to original margins.

#' Transform the columns of a dataset to (approximately) unit Frechet margins
#'
#' @description Transforms columns of dataset to unit Frechet margins, to ensure
#' the theoretical requirements are satisfied for the application of
#' \code{\link{max_stable_prcomp}} using the empirical distribution function. 
#'
#' @param data, array or vector with the data which columns are to be transformed
#' @return array or vector of same shape and type as data with the transformed data with unit
#' Frechet margins-
#' @seealso [max_stable_prcomp()], [transform_orig_margins()], [mev::fit.gev())] for information about why to transform data. 
#' @export
#' @examples
#' # sample some data
#' dat <- rnorm(1000)
#' transformed_dat <- transform_unitfrechet(dat)
#'
#' # Look at a plot of distribution
#' boxplot(transformed_dat)
#' plot(stats::ecdf(transformed_dat))
transform_unitfrechet <- function(data) {
  # turn dataframes into a matrix
  if (is.data.frame(data)) data <- as.matrix(data)

  if (!is.array(data)) {
    n <- length(data)
    d <- 1
  } else {
    n <- dim(data)[1]
    d <- dim(data)[2]
  }

  # Inspired by SpatialExtremes::gev2frech
  ecdf_col <- function(x) stats::ecdf(x)(x) * n / (n + 1)

  # check for dimensions
  if (d == 1) {
    data_ecdf <- ecdf_col(data)
    data_unitfrech <- -1 / log(data_ecdf)
  } else {
    data_ecdf <- apply(data, 2, ecdf_col)
    data_unitfrech <- apply(data_ecdf, 2, function(x) -1 / log(x))
  }

  return(data_unitfrech)
}

#' Transform the columns of a dataset to unit Pareto
#'
#' @description Transforms columns of dataset to unit Pareto margins, to ensure
#' the theoretical requirements are satisfied for the application of
#' \code{\link{max_stable_prcomp}} using the empirical distribution function. 
#'
#' @param data, array or vector with the data which columns are to be transformed
#' @return array or vector of same shape and type as data with the transformed data with unit
#' Frechet margins-
#' @seealso [max_stable_prcomp()], [transform_orig_margins()], [mev::fit.gev())] for information about why to transform data. 
#' @export
#' @examples
#' # sample some data
#' dat <- rnorm(1000)
#' transformed_dat <- transform_unitfrechet(dat)
#'
#' # Look at a plot of distribution
#' boxplot(transformed_dat)
#' plot(stats::ecdf(transformed_dat))
transform_unitpareto <- function(data) {
  # turn dataframes into a matrix
  if (is.data.frame(data)) data <- as.matrix(data)

  if (!is.array(data)) {
    n <- length(data)
    d <- 1
  } else {
    n <- dim(data)[1]
    d <- dim(data)[2]
  }

  # Inspired by de Haan Einmahl Piterbarg
  ecdf_col <- function(x) stats::ecdf(x)(x) - 1 / n

  # check for dimensions
  if (d == 1) {
    data_ecdf <- ecdf_col(data)
    data_unitpareto <- 1 / (1 - data_ecdf)
  } else {
    data_ecdf <- apply(data, 2, ecdf_col)
    data_unitpareto <- apply(data_ecdf, 2, function(x) 1 / (1 - ecdf_col(x)))
  }

  return(data_unitpareto)
}


#' Transform the columns of a transformed dataset to original margins
#'
#' Since the dataset is intended to be transformed for PCA,
#' this function takes a dataset transformed_data and
#' transforms the margins to the marginal distribution of
#' the dataset orig_data.
#'
#' @param transformed_data, arraylike data of dimension n, d
#' @param orig_data, arraylike data of dimension n , d
#' @return array of dimension n,d with transformed columns of transformed_data that follow approximately the same
#' marginal distribution of orig_data.
#' @seealso [max_stable_prcomp()], [transform_unitfrechet()], [mev::fit.gev())] for information about why to transform data
#' @export
#' @examples
#' # create a sample
#' dat <- rnorm(1000)
#' transformed_dat <- transform_unitpareto(dat)
transform_orig_margins <- function(transformed_data, orig_data) {
  # turn dataframes into a matrix
  if (is.data.frame(orig_data)) orig_data <- as.matrix(orig_data)
  if (is.data.frame(transformed_data)) transformed_data <- as.matrix(transformed_data)

  # check for dimension of data
  if (!is.array(transformed_data)) {
    n <- length(transformed_data)
    d <- 1
  } else {
    n <- dim(transformed_data)[1]
    d <- dim(transformed_data)[2]
  }

  # set up result matrix
  result <- matrix(0, n, d)

    if (d == 1) {
      qresult <- stats::quantile(orig_data, exp(-1 / transformed_data))
      result[, 1] <- as.vector(qresult)
    } else {
      for (j in 1:d) {
        # calculate transformation by empircal pseudoinverse of original data applied to uniform
        # distribution obtained by applying distribution function to standardized data.
        qresult <- stats::quantile(orig_data[, j], exp(-1 / transformed_data[, j]))

        # Quantile function returns not really a vector so transform to compatible type
        # to set matrix column.
        result[, j] <- as.vector(qresult)
      }
    }
  return(result)
}