# --- main function to calculate max_stable_prcomp object ---

#' Calculate max-stable PCA with dimension p for given dataset
#'
#' Find a p dimensional max-linear subspace that describes the data 
#' best by trying to find optimal max-linear combinations. 
#' 
#' For details on the statistical procedure it is advised to 
#' read the articles TODO. 
#'
#' @param data, array or data.frame of n observations of d variables
#' with unit Frechet margins. The max-stable PCA is fitted to
#' reconstruct this dataset with a rank p approximation.
#' @param p, integer between 1 and ncol(data). Determines
#' the dimension of the encoded state, i.e. the number of max-linear
#' combinations in the compressed representation.
#' @param n_initial_guesses number of guesses to choose a valid initial value 
#' for optimization from. This procedure uses a pseudo random number generator so 
#' setting a seed is necessary for reproducibility. 
#' @param s (default = 3), numeric greater than 0. Hyperparameter for the
#' stable tail dependence estimator used in tn the calculation.
#' @param ... additional parameters passed to \code{link{nloptr::slsqp()}}
#' @return object of class max_stable_prcomp with slots
#' p, inserted value of dimension,
#' decoder_matrix, an array of shape (d,p), where the columns represent the basis of the max-linear space for the reconstruction.
#' encoder_matrix, an array of shape (p,d), where the rows represent the loadings as max-linear combinations for the compressed representation.
#' reconstr_matrix, an array of shape (d,d), where the matrix is the mapping of the data to the optimal max-linear combinations of columns of decoder_matrix.
#' loss_fctn_value, float representing the final loss function value of the fit.
#' optim_conv_status, integer indicating the convergence of the optimizer if greater than 0.
#' @export
#' @examples
#' # generate some data with the desired margins
#' dat <- matrix(evd::rfrechet(300), 100, 3)
#' maxPCA <- max_stable_prcomp(dat, 2)
#' 
#' # look at summary to obtain further information about 
#' # loadings the space spanned and loss function
#' summary(maxPCA)
#' 
#' # transfrom data to compressed representation
#' # for a representation that is p-dimensional,
#' # preserves the max-stable structure and is numeric solution to 
#' # optimal reconstruction.
#' compr <- compress(maxPCA, dat)
#' 
#' # For visual examination reconstruct original vector from compressed representation
#' rec <- reconstruct(maxPCA, compr)
max_stable_prcomp <- function(data, p, s = 3, n_initial_guesses = 150, ...) {

  data <- as.matrix(data)
  d <- dim(data)[2]

  # setting up target function with data to be callable
  target_fn <- function(x) target_fn_data(x, d, p, s, data)

  # setting up lower bound to absolutely prohibit negative values
  lower <- rep(0, 2 * d * p)

  # setting up sharper inequality constraints than just lower bound
  constr_hin <- function(x) constraint_ineq(x, d, p)

  x0 <- NA

  # norm data and calculate distances between points from data on 1-norm sphere
  data_normed <- t(apply(data, 1, function(x) x / sum(abs(x))))
  distmat <- matrix(0, nrow(data_normed), nrow(data_normed))
  for(i in 1:nrow(data_normed)) {
    for(j in 1:i) {
      distmat[i,j] <- sum(abs(data_normed[i,] - data_normed[j,]))
    }
  }

  # search for p data entries which have a big distance from one another
  counter <- 0
  indices_base <- c()
  while(counter < p) {
    indices <- which(distmat == max(distmat), arr.ind = T)
    distmat[indices[1], indices[2]] <- 0
    for(ind in indices) {
      if(!(ind %in% indices_base) & length(indices_base) < p) {
        indices_base <- c(indices_base, ind)
        counter <- counter + 1
      }
    }
  }


  # create the base columns for the first d*p entries of inital value
  data_entries <- data[indices_base, ]
  base_columns <- apply(data_entries, 1, function(x) x / max(abs(x)))

  # search for best starting point with the specified max-space 
  # filling up the remaining d*p columns with uniform entries
  searching_x0 <- T

  while(searching_x0) {
    x0_cands <- matrix(stats::runif(n_initial_guesses * 2 * d * p, 0.1, 1.15), n_initial_guesses, 2 * d * p)
    x0_cands[, 1:(d*p)] <- matrix(base_columns, n_initial_guesses, d * p, byrow = T)

    x0_valid <- x0_cands[apply(x0_cands, 1, function(x) all(constr_hin(x) >= 0)), ]
    if(length(x0_valid) > 0) {
      searching_x0 <- F
      if(length(x0_valid) > 1) {
        targetvals <- apply(x0_valid, 1, target_fn)
        print(paste("Number of valid inits:", length(targetvals)))
        x0 <- x0_valid[which(targetvals == min(targetvals)), ]
      } else {
        x0 <- x0_valid
      }
    }
  }

# call the optimizer to calculate a solution candidate
optimizer_result <- nloptr::slsqp(x0, target_fn, hin = constr_hin, lower = lower, ...)

# set up the necessary matrices and objects for the return value
decoder_matrix <- matrix(optimizer_result$par[1:(d*p)], d, p)
encoder_matrix <- matrix(optimizer_result$par[(d*p + 1):(2 * d * p)], p, d)
reconstr_matrix <- maxmatmul(decoder_matrix, encoder_matrix)

result <- list(
               p = p,
               d = nrow(reconstr_matrix),
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
#' # generate some data with the desired margins
#' dat <- matrix(evd::rfrechet(300), 100, 3)
#' maxPCA <- max_stable_prcomp(dat, 2)
#' 
#' #  look at summary to obtain further information about 
#' # loadings the space spanned and loss function
#' summary(maxPCA)
#' 
#' # transfrom data to compressed representation
#' # for a representation that is p-dimensional,
#' # preserves the max-stable structure and is numeric solution to 
#' # optimal reconstruction.
#' compr <- compress(maxPCA, dat)
#' 
#' # For visual examination reconstruct original vector from compressed representation
#' rec <- reconstruct(maxPCA, compr)
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
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @param compressed_data a compressed representation of 
#' data obtained by calling \code{compress(fit, data)} 
#' intended to be calculated with the same fit.
#' @return array of dimension nrow(compressed_data), d
#' An optimal reconstruction with max-linear combinations 
#' from p max-linearly independent columns.
#' @seealso [maxstablePCA::compress(), 
#' maxstablePCA::max_stable_prcomp(),
#' maxstablePCA::maxmatmul()]
#' @export
#' @examples
#' # generate some data with the desired margins
#' dat <- matrix(evd::rfrechet(300), 100, 3)
#' maxPCA <- max_stable_prcomp(dat, 2)
#' 
#' #  look at summary to obtain further information about 
#' # loadings the space spanned and loss function
#' summary(maxPCA)
#' 
#' # transfrom data to compressed representation
#' # for a representation that is p-dimensional,
#' # preserves the max-stable structure and is numeric solution to 
#' # optimal reconstruction.
#' compr <- compress(maxPCA, dat)
#' 
#' # For visual examination reconstruct original vector from compressed representation
#' rec <- reconstruct(maxPCA, compr)
reconstruct <- function(fit, compressed_data) {
  # TODO: take care of vectors!
  return(t(maxmatmul(fit$decoder_matrix, t(compressed_data))))
}

#' Print summary of a max_stable_prcomp object.
#'
#' 
#'
#' @param object, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @param ... additional unused arguments.
#' @export
#' @examples
#' # Alternatively call via the generic function summary
#' # TODO: generic fit example
#' # summary.max_stable_prcomp(fit)
summary.max_stable_prcomp <- function(object, ...) {
  print(object)

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
  # remove lower constraint for the random initialization
  # constr_lower <- rowSums(H) - 1
  return(constr_upper)
}

