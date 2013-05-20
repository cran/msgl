\name{msgl.cv}
\alias{msgl.cv}
\title{Multinomial sparse group lasso cross validation using multiple possessors}
\usage{
  msgl.cv(x, classes, sampleWeights = NULL,
    grouping = NULL, groupWeights = NULL,
    parameterWeights = NULL, alpha = 0.5,
    standardize = TRUE, lambda, fold = 10L,
    cv.indices = list(), sparse.data = FALSE,
    max.threads = 2L, seed = 331L,
    algorithm.config = sgl.standard.config)
}
\arguments{
  \item{x}{design matrix, matrix of size \eqn{N \times p}.}

  \item{classes}{classes, factor of length \eqn{N}.}

  \item{sampleWeights}{sample weights, a vector of length
  \eqn{N}.}

  \item{grouping}{grouping of covariates, a vector of
  length \eqn{p}. Each element of the vector specifying the
  group of the covariate.}

  \item{groupWeights}{the group weights, a vector of length
  \eqn{m+1} (the number of groups). The first element of
  the vector is the intercept weight. If \code{groupWeights
  = NULL} default weights will be used. Default weights are
  0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number
  of covariates in the group}}} for all other weights.}

  \item{parameterWeights}{a matrix of size \eqn{K \times
  (p+1)}. The first column of the matrix is the intercept
  weights. Default weights are is 0 for the intercept
  weights and 1 for all other weights.}

  \item{alpha}{the \eqn{\alpha} value 0 for group lasso, 1
  for lasso, between 0 and 1 gives a sparse group lasso
  penalty.}

  \item{standardize}{if TRUE the covariates are standardize
  before fitting the model. The model parameters are
  returned in the original scale.}

  \item{lambda}{the lambda sequence for the regularization
  path.}

  \item{fold}{the fold of the cross validation, an integer
  larger than \eqn{1} and less than \eqn{N+1}. Ignored if
  \code{cv.indices != NULL}. If
  \code{fold}\eqn{\le}\code{max(table(classes))} then the
  data will be split into \code{fold} disjoint subsets
  keeping the ration of classes approximately equal.
  Otherwise the data will be split into \code{fold}
  disjoint subsets without keeping the ration fixed.}

  \item{cv.indices}{a list of indices of a cross validation
  splitting. If \code{cv.indices = NULL} then a random
  splitting will be generated using the \code{fold}
  argument.}

  \item{sparse.data}{if TRUE \code{x} will be treated as
  sparse, if \code{x} is a sparse matrix it will be treated
  as sparse by default.}

  \item{max.threads}{the maximal number of threads to be
  used}

  \item{seed}{the seed used for generating the random cross
  validation splitting, only used if
  \code{fold}\eqn{\le}\code{max(table(classes))}.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{link}{the linear predictors -- a list of length
  \code{length(lambda)} one item for each lambda value,
  with each item a matrix of size \eqn{K \times N}
  containing the linear predictors.} \item{response}{the
  estimated probabilities - a list of length
  \code{length(lambda)} one item for each lambda value,
  with each item a matrix of size \eqn{K \times N}
  containing the probabilities.} \item{classes}{the
  estimated classes - a matrix of size \eqn{N \times d}
  with \eqn{d=}\code{length(lambda)}.}
  \item{cv.indices}{the cross validation splitting used.}
  \item{features}{average number of features used in the
  models.} \item{parameters}{average number of parameters
  used in the models.}
}
\description{
  Multinomial sparse group lasso cross validation using
  multiple possessors
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.03)
fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)

# Missclassification count
colSums(fit.cv$classes != classes)
}
\author{
  Martin Vincent
}

