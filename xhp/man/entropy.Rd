% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy.R
\name{entropy}
\alias{entropy}
\title{Clustering Entropy}
\usage{
entropy(x)
}
\arguments{
\item{x}{A vector of cluster labels}
}
\value{
A scalar giving the clustering entropy.
}
\description{
The clustering entropy is computed.  Let \eqn{p_i} be the relative size of
cluster \eqn{i}. The clustering entropy is the negative of the sum of
\eqn{p_i * f(p_i)} for all \eqn{i}, where \eqn{f} is the natural logarithm.
}
\examples{
entropy(c(1, 1, 1, 1, 1))
entropy(c(1, 1, 2, 2, 2))
entropy(c(1, 2, 3, 4, 5))

}
