% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{transform_orig_margins}
\alias{transform_orig_margins}
\title{Transform the columns of a transformed dataset to original margins}
\usage{
transform_orig_margins(transformed_data, orig_data)
}
\arguments{
\item{transformed_data, }{arraylike data of dimension n, d}

\item{orig_data, }{arraylike data of dimension n , d}
}
\value{
array of dimension n,d with transformed columns of transformed_data that follow approximately the same
marginal distribution of orig_data.
}
\description{
Since the dataset is intended to be transformed for PCA,
this function takes a dataset transformed_data and
transforms the margins to the marginal distribution of
the dataset orig_data.
}
\examples{
# create a sample
dat <- rnorm(1000)
transformed_dat <- transform_unitpareto(dat)
}
\seealso{
\code{\link[=max_stable_prcomp]{max_stable_prcomp()}}, \code{\link[=transform_unitfrechet]{transform_unitfrechet()}}, [mev::fit.gev())] for information about why to transform data

[mev::fit.gev())]: R:mev::fit.gev())
}
