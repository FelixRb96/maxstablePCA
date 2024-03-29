% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{reconstruct}
\alias{reconstruct}
\title{Reconstruct the data from the compressed representation}
\usage{
reconstruct(fit, compressed_data)
}
\arguments{
\item{fit, }{max_stable_prcomp object. Data should be
assumed to follow the same distribution as the data used in
max_stable_prcomp.}

\item{compressed_data}{a compressed representation of 
data obtained by calling \code{compress(fit, data)} 
intended to be calculated with the same fit.}
}
\value{
array of dimension nrow(compressed_data), d
An optimal reconstruction with max-linear combinations 
from p max-linearly independent columns.
}
\description{
Calculate the reconstruction given by the compressed state
compressed_data and the fit object from the max-stable PCA.
Can be used for visual assessment of the fit if p should be chosen
differently.
}
\examples{
# generate some data with the desired margins
dat <- matrix(evd::rfrechet(300), 100, 3)
maxPCA <- max_stable_prcomp(dat, 2)

#  look at summary to obtain further information about 
# loadings the space spanned and loss function
summary(maxPCA)

# transfrom data to compressed representation
# for a representation that is p-dimensional,
# preserves the max-stable structure and is numeric solution to 
# optimal reconstruction.
compr <- compress(maxPCA, dat)

# For visual examination reconstruct original vector from compressed representation
rec <- reconstruct(maxPCA, compr)
}
\seealso{
[maxstablePCA::compress(), 
maxstablePCA::max_stable_prcomp(),
maxstablePCA::maxmatmul()]
}
