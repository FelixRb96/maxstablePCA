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
#' @param norm (delfault "l1") which norm to use for the spectral measure estimator, currently only l1 and sup norm "linfty" are available. 
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
max_stable_prcomp <- function(data, p, s = 3, lambda = 0, norm = "l1", n_initial_guesses = 150, ...) {

  dat <- as.matrix(data)
  if(norm == "l1") {
    dat_extr <- data[which(rowSums(dat) > s), ]
    dat_normed_extr <- t(apply(dat_extr, 1, function(z) z / sum(abs(z))))
  }
  if(norm == "linfty") {
    dat_extr <- data[which(apply(dat, 1, max) > s), ]
    dat_normed_extr <- t(apply(dat_extr, 1, function(z) z / max(z)))
  }



  n <- dim(data)[1]
  d <- dim(data)[2]
  
  target_fn <- function(x) target_fn_data(x, d, p, s, n, lambda, dat_normed_extr)

  # setting up lower bound to absolutely prohibit negative values
  lower <- rep(0, 2 * d * p)
  upper <- c(rep(as.numeric(d - p + 1),d * p), rep(1.0, d * p))


  # setting up sharper inequality constraints than just lower bound
  constr_hin <- function(x) constraint_ineq(x, d, p)

  x0 <- NA


  # take care of the case p = 1 later where distances are not relevant
  if(p > 1) {

    # make sure the distance matrix does not get too large
    if(nrow(dat_extr) > 5000) {
      warning("s seems to be small or the data set is large, optimization might fail...")
    }


    # norm data and calculate distances between points from data on 1-norm sphere
    distmat <- matrix(0, nrow(dat_normed_extr), nrow(dat_normed_extr))
    for(i in 1:nrow(dat_normed_extr)) {
      for(j in 1:i) {
        distmat[i,j] <- sum(abs(dat_normed_extr[i,] - dat_normed_extr[j,]))
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
      x0_cands <- matrix(stats::runif(n_initial_guesses * 2 * d * p, 0.1, 1), n_initial_guesses, 2 * d * p)
      x0_cands[, 1:(d*p)] <- matrix(base_columns, n_initial_guesses, d * p, byrow = T)

      x0_valid <- x0_cands[apply(x0_cands, 1, function(x) all(constr_hin(x) >= 0)), ]
      if(length(x0_valid) > 0) {
        searching_x0 <- F
        if(length(x0_valid) > 1) {
          targetvals <- apply(x0_valid, 1, target_fn)
          x0 <- x0_valid[which(targetvals == min(targetvals)), ]
        } else {
          x0 <- x0_valid
        }
      }
    }
  } else {
    x0 <- c(base::colMeans(dat_extr) / sum(base::colMeans(dat_extr)), rep(0.999, d))
  }

  # call the optimizer to calculate a solution candidate
  # opts <- list("algorithm"="NLOPT_LN_COBYLA")
  # optimizer_result <- nloptr::nloptr(x0, target_fn, eval_g_ineq = function(x) -constr_hin(x), lb = lower, opts = opts, ...)
  optimizer_result <- nloptr::slsqp(
                                    x0, 
                                    target_fn, 
                                    hin = constr_hin, 
                                    lower = lower, 
                                    upper = upper, 
                                    ...
  )

  # set up the necessary matrices and objects for the return value
  # decoder_matrix <- matrix(optimizer_result$solution[1:(d*p)], d, p)
  # encoder_matrix <- matrix(optimizer_result$solution[(d*p + 1):(2 * d * p)], p, d)
  encoder_matrix <- matrix(optimizer_result$par[(d*p + 1):(2 * d * p)], p, d)
  decoder_matrix <- matrix(optimizer_result$par[1:(d*p)], d, p)

  reconstr_matrix <- maxmatmul(decoder_matrix, encoder_matrix)

  result <- list(
                 p = p,
                 d = nrow(reconstr_matrix),
                 decoder_matrix = decoder_matrix,
                 encoder_matrix = encoder_matrix,
                 reconstr_matrix = reconstr_matrix,
                 loss_fctn_value = optimizer_result$value,
                 optim_conv_status = optimizer_result$convergence,
                 s = s, 
                 starting_vals = list(
                                      encoder_matrix_x0 = matrix(x0[1:(d * p)], d, p), 
                                      decoder_matrix_xo = matrix(x0[(d * p + 1):(2 * d * p)], d, p)
                 )
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

#' Loss for a given reconstruction matrix
#' 
#' TODO 
#' 
#' @param H, a (d x d) matrix where d is the number of columns of data
#' @param data, a matrix like object of shape (n x d), where d is the same as for H. Should be transformed (TODO). 
#' @param s, non-negative numeric, tunig parameter for the stable tail dependence function. For further information see TODO. 
#' @export
#' @examples 
#' TODO
reconstr_loss <- function(H, data, s) {

  id <- diag(length(H[1,]))

  std <- function(z) stable_tail_dependenceC(z, s, data)
  row_fctn <- function(i) return(2 * std(pmax(H[i,], id[i,])) - std(H[i,]))

  # calculate result
  result <- sum(sapply(1:length(H[1,]), row_fctn))
  return(result)

}



# --- internal functions to solve minimization problem ---

# Creates the reconstruction matrix H from optimizer result
create_H <- function(x, d, p) {
  A <- matrix(x[1:(d * p)], d, p)
  V <- matrix(x[(d * p + 1):length(x)], p, d)
  H <- maxmatmul(A, V)
  return(H)
}

target_fn_data <- function(x, d, p, s, n, lambda, data_normed_extr) {
  # create H from optimizer vector
  H <- create_H(x, d, p)

  # calculate rownorms and change to relevant subset
  rec_standardized <- t(maxmatmul(H, t(data_normed_extr)))
  penalty <- lambda * (sum(abs(x)))
  return(sum(abs(data_normed_extr - rec_standardized)) * s / n + lambda * penalty)
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

