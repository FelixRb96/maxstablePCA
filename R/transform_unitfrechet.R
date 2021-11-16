# Tools for transforming data to the necessary margins
# and given an analyzed transformed dataset the possiblity
# to transform to original margins.

#' Transform the columns of a dataset to unit Frechet margins
#'
#' Lore ipsum
#'
#' @param data, array or vector with the data which columns are to be transformed
#' @param empirical, logical to indicate the transformation method. If true the data
#' is transformed using the empical distribution function. If \code{FALSE}, the 
#' data is transformed using a MLE fit of a generalized extreme value distribution to the data.
#' Using the empirical distribution function does not require the data to have approximately
#' max-stable distribution, thus is set to default. 
#' @return array or vector of same shape and type as data with the transformed data with unit 
#' Freceht margins-
#' @seealso [maxstablePCA::max_stable_prcomp()] for information about why to transform data.
#' @export
#' @examples
#' # TODO
#' 1 + 1
transform_unitfrechet <- function(data, empirical = TRUE) {
    if(is.vector(data)) {
      n <- length(data)
      d <- 1
    } else {
      n <- dim(data)[1]
      d <- dim(data)[2]
    }

  if(empirical) {

    # Inspired by SpatialExtremes::gev2frech (Thanks Mr. Ribatet!)
    ecdf_col <- function(x) ecdf(x)(x) * n / (n + 1)

    # check for dimensions
    if(d == 1) {
      data_ecdf <- ecdf_col(data)
      data_unitfrech <- -1 / log(data_ecdf)
    } else {
        data_ecdf <- apply(data, 2, ecdf_col)
        data_unitfrech <- apply(data_ecdf, 2, function(x) -1 / log(x))
    }

    return(data_unitfrech)
  } else {

    # transformation function
    transform_col_frech <- function(x, loc, scale, shape) return(pmax(1 + shape * ( ( x - loc) / scale),0)^(1 / shape))

    # check for dimensions
    if(d == 1) {
      fits <- mev::fit.gev(data)
      data_unitfrech <- transform_col_frech(data, fits$estimate[1], fits$estimate[2], fits$estimate[3])
    } else {
      fits <- apply(data, 2, mev::fit.gev)
      transform_by_index <- function(i) transform_col_frech(data[,i], fits[[i]]$estimate[1], fits[[i]]$estimate[2], fits[[i]]$estimate[3])
      data_unitfrech <- sapply(1:d, transform_by_index)
    }

    return(data_unitfrech)
  }
}

#' Transform the columns of a transformed dataset to original margins
#'
#' Lore ipsum
#'
#' @param A TODO
#' @param B TODO
#' @return TODO
#' @seealso [maxstablePCA::max_stable_prcomp()] for information about why to transform data.
#' @export
#' @examples
#' # TODO
#' 1 + 1
transform_orig_margins <- function(transformed_data, orig_data, empirical = TRUE) {
  return("TODO")
}


