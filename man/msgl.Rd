\name{msgl}
\alias{msgl}
\title{Fit a multinomial sparse group lasso regularization path.}
\usage{
  msgl(x, classes,
    sampleWeights = rep(1/length(classes), length(classes)),
    grouping = NULL, groupWeights = NULL,
    parameterWeights = NULL, alpha = 0.5,
    standardize = TRUE, lambda, return = 1:length(lambda),
    sparse.data = FALSE,
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

  \item{return}{the indices of lambda values for which to
  return a the fitted parameters.}

  \item{sparse.data}{if TRUE \code{x} will be treated as
  sparse, if \code{x} is a sparse matrix it will be treated
  as sparse by default.}

  \item{algorithm.config}{the algorithm configuration to be
  used.}
}
\value{
  \item{beta}{the fitted parameters -- a list of length
  \code{length(lambda)} with each entry a matrix of size
  \eqn{K\times (p+1)} holding the fitted parameters}
  \item{loss}{the values of the loss function}
  \item{objective}{the values of the objective function
  (i.e. loss + penalty)} \item{lambda}{the lambda values
  used}
}
\description{
  For a classification problem with \eqn{K} classes and
  \eqn{p} covariates dived into \eqn{m} groups. A sequence
  of minimizers (one for each lambda given in the
  \code{lambda} argument) of \deqn{\hat R(\beta) + \lambda
  \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2
  + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)} where
  \eqn{\hat R} is the weighted empirical log-likelihood
  risk of the multinomial regression model. The vector
  \eqn{\beta^{(J)}} denotes the parameters associated with
  the \eqn{J}'th group of covariates (default is one
  covariate per group, hence the default dimension of
  \eqn{\beta^{(J)}} is \eqn{K}). The group weights
  \eqn{\gamma \in [0,\infty)^m} and the parameter weights
  \eqn{\xi = (\xi^{(1)},\dots, \xi^{(m)}) \in [0,\infty)^n}
  with \eqn{\xi^{(1)}\in [0,\infty)^{n_1},\dots, \xi^{(m)}
  \in [0,\infty)^{n_m}}.
}
\examples{
data(SimData)
x <- sim.data$x
classes <- sim.data$classes
lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
fit <- msgl(x, classes, lambda = lambda)
fit$beta[[10]] #model with lambda = lambda[10]
}
\author{
  Martin Vincent
}

