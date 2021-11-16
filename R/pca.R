# --- main function to calculate max_stable_prcomp object ---

#' Calculate max-stable PCA with dimension p for given dataset
#'
#' Find
#'
#' @param data, array or data.frame of n observations of d variables
#' with unit Frechet margins. The max-stable PCA is fitted to
#' reconstruct this dataset with a rank p approximation.
#' @param p, integer between 1 and ncol(data). Determines
#' the dimension of the encoded state, i.e. the number of max-linear
#' combinations in the compressed representation.
#' @param s (default = 3), numeric greater than 0. Hyperparameter for the
#' stable tail dependence estimator used in tn the calculation.
#' @param x0 (default = NULL), initial value for optimizer. Set with
#' extreme caution and only with good knowledge of the minimizationb problem.
#' NOTE: will be moved to ...
#' @return object of class max_stable_prcomp with slots
#' TODO
#' @export
#' @examples
#' # TODO
#' 1 + 1
max_stable_prcomp <- function(data, p, s = 3, x0 = NULL) {

  data <- as.matrix(data)
  d <- dim(data)[2]

  # setting up an intial guess thats always feasible
  if(is.null(x0)) {
    x0 <- rep(0, 2 * d * p)
    for(i in 1:p) {
      x0[(i-1)*d + i] <- 1
      x0[(d * (i-1) + (p+1)):(d*i)] <- 1/p
      x0[d*p + (i-1)*p + i] <- 1
    }
  }

  # setting up target function with data to be callable
  target_fn <- function(x) target_fn_data(x, d, p, s, data)

  # setting up lower bound to absolutely prohibit negative values
  lower <- rep(0, 2 * d * p)

  # setting up sharper inequlity constraints
  constr_hin <- function(x) constraint_ineq(x, d, p)


  optimizer_result <- nloptr::slsqp(x0, target_fn, hin = constr_hin, lower = lower)

  decoder_matrix <- matrix(optimizer_result$par[1:(d*p)], d, p)
  encoder_matrix <- matrix(optimizer_result$par[(d*p + 1):(2 * d * p)], p, d)
  reconstr_matrix <- maxmatmul(decoder_matrix, encoder_matrix)

  result <- list(
                 p = p,
                 decoder_matrix = decoder_matrix,
                 encoder_matrix = encoder_matrix,
                 reconstr_matrix = reconstr_matrix,
                 loss_fctn_value = optimizer_result$value - d,
                 optim_conv_status = optimizer_result$convergence,
                 s = s
  )

  class(result) <- "max_stable_prcomp"
  return(result)
}

# --- helper functions to obtain latent space representations and transofrm data ---

#' Transform data to compact representation given by max-stable PCA
#'
#' Turn the given data into a compressed latent representation
#' given by the fit of the max_stable_prcomp function.
#' This is done by taking the max-matrix product of the data
#' and the encoder matrix from the fit.
#'
#' @param data, array with same number of columns as the
#' data of the fit object.
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @seealso [maxstablePCA::max_stable_prcomp(), maxstablePCA::maxmatmul()]
#' @return An array of shape nrow(data), p giving the
#' encoded representation of the data in p components which are
#' also unit Frechet distributed which is to be takin into consideration for
#' further analysis.
#' @export
#' @examples
#' # TODO
#' 1 + 1
compress <- function(fit, data) {
  # TODO: take care of vectors!
  return(t(maxmatmul(fit$encoder_matrix, t(data))))
}

#' Reconstruct the data from the compressed representation
#'
#' Calculate the reconstruction given by the compressed state
#' compressed_data and the fit object from the max-stable PCA.
#' Can be used for visual assessment of the fit if p should be chosen
#' differently.
#'
#'
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @param compressed_data a compressed representation of 
#' data obtained by calling \code[compress(fit, data)} 
#' intended to be calculated with the same fit.
#' @return array of dimension nrow(compressed_data), d
#' An optimal reconstruction with max-linear combinations 
#' from p max-linearly independent columns.
#' @seealso [maxstablePCA::compress(), 
#' maxstablePCA::max_stable_prcomp(),
#' maxstablePCA::maxmatmul()]
#' @export
#' @examples
#' # TODO
#' 1 + 1
reconstruct <- function(fit, compressed_data) {
  # TODO: take care of vectors!
  return(t(maxmatmul(fit$decoder_matrix, t(compressed_data))))
}

#' Print summary of a max_stable_prcomp object.
#'
#' TODO
#'
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @export
#' @examples
#' # Alternatively call via the generic function summary
#' # TODO: generic fit example
#' # summary.max_stable_prcomp(fit)
summary.max_stable_prcomp <- function(fit) {
  # TODO: nice formatting with cat or sth.
  print(fit)
}



# --- internal functions to solve minimization problem ---

# Creates the reconstruction matrix H from optimizer result
create_H <- function(x, d, p) {
  A <- matrix(x[1:(d * p)], d, p)
  V <- matrix(x[(d * p + 1):length(x)], p, d)
  H <- maxmatmul(A, V)
  return(H)
}

# Helper function to calculate the distance of the reconstruction 
#  from the data.
target_fn_data <- function(x, d, p, s, data)  {

  # set up matrices and functions needed from input
  H <- create_H(x, d, p)
  id <- diag(length(H[1,]))

  std <- function(z) stable_tail_dependence(z, s, data)
  row_fctn <- function(i) return(2 * std(pmax(H[i,], id[i,])) - std(H[i,]))

  # calculate result
  result <- sum(sapply(1:length(H[1,]), row_fctn))
  return(result)
}

# constraint helper function to ensure non-neagtivity of vector x.
constraint_nonneg <- function(x) return(x)

# Constraint helper to ensure the bounds of the rows of H.
constraint_ineq <- function(x, d, p) {
  H <- create_H(x, d, p)
  constr_upper <- d - rowSums(H)
  constr_lower <- rowSums(H) - 1
  return(c(constr_lower, constr_upper))
}

