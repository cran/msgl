\name{msgl.subsampling}
\alias{msgl.subsampling}
\title{Multinomial sparse group lasso generic subsampling procedure}
\usage{
  msgl.subsampling(x, classes,
    sampleWeights = rep(1/length(classes), length(classes)),
    grouping = NULL, groupWeights = NULL,
    parameterWeights = NULL, alpha = 0.5,
    standardize = TRUE, lambda, training, test,
    sparse.data = FALSE, max.threads = 2L,
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

  \item{training}{a list of training samples, each item of
  the list corresponding to a subsample. Each item in the
  list must be a vector with the indices of the training
  samples for the corresponding subsample. The length of
  the list must equal the length of the \code{test} list.}

  \item{test}{a list of test samples, each item of the list
  corresponding to a subsample. Each item in the list must
  be vector with the indices of the test samples for the
  corresponding subsample. The length of the list must
  equal the length of the \code{training} list.}

  \item{sparse.data}{if TRUE \code{x} will be treated as
  sparse, if \code{x} is a sparse matrix it will be treated
  as sparse by default.}

  \item{max.threads}{the maximal number of threads to be
  used}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{link}{the linear predictors -- a list of length
  \code{length(test)} with each element of the list another
  list of length \code{length(lambda)} one item for each
  lambda value, with each item a matrix of size \eqn{K
  \times N} containing the linear predictors.}
  \item{response}{the estimated probabilities -- a list of
  length \code{length(test)} with each element of the list
  another list of length \code{length(lambda)} one item for
  each lambda value, with each item a matrix of size \eqn{K
  \times N} containing the probabilities.}
  \item{classes}{the estimated classes -- a list of length
  \code{length(test)} with each element of the list a
  matrix of size \eqn{N \times d} with
  \eqn{d=}\code{length(lambda)}.} \item{features}{number of
  features used in the models.} \item{parameters}{number of
  parameters used in the models.}
}
\description{
  Support the use of multiple processors.
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.03)

test <- replicate(5, sample(1:length(classes))[1:20], simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda,
 training = train, test = test)

# Missclassification count of second subsample
colSums(fit.sub$classes[[2]] != classes[test[[2]]])
}
\author{
  Martin Vincent
}

