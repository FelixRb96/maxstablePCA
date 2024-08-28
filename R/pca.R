# --- main function to calculate max_stable_prcomp object ---

#' Calculate max-stable PCA with dimension p for given dataset
#'
#' @description Find an optimal encoding of data of extremes using max-linear combinations
#' by a distance minimization approach. Can be used to check if the data
#' follows approximately a generalized max-linear model.
#' For details on the statistical procedure it is advised to
#' consult the articles "F. Reinbott, A. Janßen, Principal component analysis for max-stable distributions (https://arxiv.org/abs/2408.10650)"
#' and "M.Schlather F. Reinbott, A semi-group approach to Principal Component Analysis (https://arxiv.org/abs/2112.04026)". 
#'
#' @name max_stable_prcomp
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
#' reconstr_matrix, an array of shape (d,d), where the matrix is the mapping of the data to the reconstruction used for the distance minimization.
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
max_stable_prcomp <- function(data, p, s = 3, n_initial_guesses = 150, norm = "l1", ...) {

  dat <- as.matrix(data)
  
  if(any(dat < 0)) base::stop("ERROR: please make sure the data is transformed to pareto margins, some values the data were negative.")

  # filter out the extreme data and prepare the normed data of extremes
  if(norm == "l1") {
    dat_extr <- data[which(rowSums(dat) > s), ]
    dat_normed_extr <- t(apply(dat_extr, 1, function(z) z / sum(abs(z))))
  } else if( norm == "linfty") {
    dat_extr <- data[which(apply(dat, 1, max) > s), ]
    dat_normed_extr <- t(apply(dat_extr, 1, function(z) z / max(abs(z))))
  } else {
    message("No valid norm specified")
    return("Failed")
  }

  # determine dimensionality
  n <- dim(data)[1]
  d <- dim(data)[2]

  target_fn <- function(x) target_fn_data(x, d, p, s, n, dat_normed_extr)
  dtarget_fn <- function(x) dtarget_fn_data(x, d, p, s, n, dat_normed_extr)

  # setting up sharper inequality constraints than just lower bound
  constr_ineq <- function(x) constr_fn(x, d, p)
  
  # random search for reasonable starting parameter for slsqp
  x0 <- NA
  searching_x0 <- T

  while(searching_x0) {
      x0_cands <- matrix(stats::runif(n_initial_guesses * 2 * d * p, 0.1, 1), n_initial_guesses, 2 * d * p)
      x0_cands[, 1:(d*p)] <- create_random_Bvals(n_initial_guesses, d, p)

      x0_valid <- x0_cands[apply(x0_cands, 1, function(x) all(constr_ineq(x) >= 0)), ]
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

  # run slsqp
  optimizer_result <- nloptr::slsqp(
    x0, 
    target_fn, 
    lower = rep(0, 2 * d * p), 
    upper = c(rep(d - p + 1, d * p), rep(1, d * p)), 
    ...
  )

  # set up the necessary matrices and objects for the return value
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
                                      decoder_matrix_x0 = matrix(x0[1:(d * p)], d, p), 
                                      encoder_matrix_x0 = matrix(x0[(d * p + 1):(2 * d * p)], p, d)
                 )
  )

  class(result) <- "max_stable_prcomp"
  return(result)
}


# --- helper functions to obtain latent space representations and transofrm data ---

#' Transform data to compact representation given by max-stable PCA
#'
#' @description Turn the given data into a compressed latent representation
#' given by the fit of the max_stable_prcomp function.
#' This is done by taking the max-matrix product of the data
#' and the encoder matrix from the fit.
#'
#' @name compress
#' @param data, array with same number of columns as the
#' data of the fit object.
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @seealso [max_stable_prcomp()], [maxmatmul()]
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
#' rec <- reconstruct(maxPCA, dat)
compress <- function(fit, data) {
  return(t(maxmatmul(fit$encoder_matrix, t(data))))
}

#' Obtain reconstructed data for PCA
#'
#' @description Map the data to the reconstruction 
#' given by the fit of the max_stable_prcomp function.
#' This is done by taking the max-matrix product of the data
#' and the reconstruction matrix from the fit.
#'
#' @name reconstruct
#' @param data, array with same number of columns as the
#' data of the fit object.
#' @param fit, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @seealso [max_stable_prcomp()], [maxmatmul()]
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
reconstruct <- function(fit, data) {
  return(t(maxmatmul(fit$reconstr_matrix, t(data))))
}

#' Print summary of a max_stable_prcomp object.
#'
#' @name summary.max_stable_prcomp
#' @param object, max_stable_prcomp object. Data should be
#' assumed to follow the same distribution as the data used in
#' max_stable_prcomp.
#' @param ... additional unused arguments.
#' @return Same as [base::print()].
#' @seealso [max_stable_prcomp()]
#' @export
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

target_fn_data <- function(x, d, p, s, n, data_normed_extr) {
  # create H from optimizer vector
  H <- create_H(x, d, p)

  # calculate rownorms and change to relevant subset
  rec_standardized <- t(maxmatmul(H, t(data_normed_extr)))
  return(sum(abs(data_normed_extr - rec_standardized)) * s / n )
}

dtarget_fn_data <- function(x, d, p, s, n, data_normed_extr) {
  A <- matrix(x[1:(d * p)], d, p)
  V <- matrix(x[(d * p + 1):length(x)], p, d)
  H <- maxmatmul(A, V)


  # set up some intermediate calculations 
  result <- rep(0, 2* d * p)
  rec <- maxmatmul(H, t(data_normed_extr))
  enc <- maxmatmul(V, t(data_normed_extr))

  for(j in 1:p) {
    for(k in 1:d) {
      tmp <- A[k,j] * enc[j,]
      if(any(tmp == rec[k,])) {
       result[(j - 1) * d + k] <- sum(sign(rec[k,]) * as.numeric((A[k,j] * enc[j,]) == (rec[k,])) * enc[j,]) 
      } else {
        result[(j - 1) * d + k] <- 0
      }
    }
  }

  for(l in 1:d) {
    for(j in 1:p) {
      tmpres <- 0
      for(i in 1:NROW(data_normed_extr)) {
        for(k in 1:d) {
          tmpres <- tmpres + sign(rec[k,i] - data_normed_extr[i,k]) * as.numeric(A[k,j] * V[j,l] * data_normed_extr[i,l] == rec[k,i]) * A[k,j] * data_normed_extr[i,l]
        }
      }
      result[d*p + (l - 1) * p + j] <- tmpres
    }
  }
  return(result * s / n)
}

create_random_Bvals <- function(n, d, p) {
   resultmat <- matrix(NA, n , d*p)
for(i in 1:n) {
  B <- matrix(stats::runif(d * p, 0.1, 1), d, p)
  B <- apply(B, 2, function(z) z / max(abs(z)))
  resultmat[i,] <- as.vector(B)
}
return(resultmat)
}

constr_fn <- function(x, d, p) {
  return(x)
}