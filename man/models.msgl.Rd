% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/navigate.R
\name{models.msgl}
\alias{models.msgl}
\title{Extract the fitted models}
\usage{
\method{models}{msgl}(object, index = 1:nmod(object), ...)
}
\arguments{
\item{object}{a msgl object}

\item{index}{indices of the models to be returned}

\item{...}{ignored}
}
\value{
a list of \eqn{\beta} matrices.
}
\description{
Returns the fitted models, that is the estimated \eqn{\beta} matrices.
}
\author{
Martin Vincent
}
