#
#     Description of this R script:
#     R interface for multinomial sparse group lasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2013 Martin Vincent
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @title Compute error rates
#'
#' @description
#' Compute error rates.
#' If \code{type = "rate"} then the misclassification rates will be computed.
#' If \code{type = "count"} then the misclassification counts will be computed.
#' If \code{type = "loglike"} then the negative log likelihood error will be computed.
#'
#' @param object a msgl object
#' @param data a matrix of
#' @param response a vector of classes
#' @param classes a vector of classes
#' @param type type of error rate
#' @param ... ignored
#' @return a vector of error rates
#'
#' @author Martin Vincent
#' @examples
#' data(SimData)
#'
#' x.all <- x
#' x.1 <- x[1:50,]
#' x.2 <- x[51:100,]
#' classes.all <- classes
#' classes.1 <- classes[1:50]
#' classes.2 <- classes[51:100]
#'
#' #### Fit models using x.1
#' lambda <- msgl::lambda(x.1, classes.1, alpha = .5, d = 25, lambda.min = 0.075)
#' fit <- msgl::fit(x.1, classes.1, alpha = .5, lambda = lambda)
#'
#' #### Training errors:
#'
#' # Misclassification rate
#' Err(fit, x.1)
#'
#' # Misclassification count
#' Err(fit, x.1, type = "count")
#'
#' # Negative log likelihood error
#' Err(fit, x.1, type="loglike")
#'
#' # Misclassification rate of x.2
#' Err(fit, x.2, classes.2)
#'
#' #### Do cross validation
#' fit.cv <- msgl::cv(x.all, classes.all, alpha = .5, lambda = lambda)
#'
#' #### Cross validation errors (estimated expected generalization error)
#'
#' # Misclassification rate
#' Err(fit.cv)
#'
#' # Negative log likelihood error
#' Err(fit.cv, type="loglike")
#'
#' #### Do subsampling
#' test <- list(1:20, 21:40)
#' train <- lapply(test, function(s) (1:length(classes.all))[-s])
#'
#' fit.sub <- msgl::subsampling(x.all, classes.all, alpha = .5,
#'  lambda = lambda, training = train, test = test)
#'
#' # Mean misclassification error of the tests
#' Err(fit.sub)
#'
#' # Negative log likelihood error
#' Err(fit.sub, type="loglike")
#'
#' @importFrom stats predict
#' @importFrom sglOptim Err
#' @importFrom sglOptim compute_error
#' @export
Err.msgl <- function(object, data = NULL, response = object$classes.true, classes = response, type = "rate", ... ) {

	if(is.null(classes)) stop("classes or response must be specified")

	loss <- switch(type,

		rate = list(function(x,y) mean(sapply(1:length(x), function(i) x[[i]] != y[[i]])), "classes", FALSE),
		count = list(function(x,y) sum(sapply(1:length(x), function(i) x[[i]] != y[[i]])), "classes", FALSE),

		loglike = list(function(x,y) -mean(log(sapply(1:length(x), function(i) y[[i]][x[[i]]]))), "response", TRUE),

		stop("Unknown type")
	)

	true_response <- classes

	if( ! is.null(data) ) {
		object <- predict(object, data)
	}

	return( compute_error(
		x = object,
		response_name = loss[[2]],
		true_response = true_response,
		loss = loss[[1]],
		transposed_response = loss[[3]])
	)
}

#' @title Nonzero features
#'
#' @description
#' Extracts the nonzero features for each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a list of of length \code{nmod(x)} containing the nonzero features (that is nonzero colums of the beta matrices)
#'
#' @examples
#' data(SimData)
#'
#'
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl::fit(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero features of model 1, 10 and 25
#' features(fit)[c(1,10,25)]
#'
#' # count the number of nonzero features in each model
#' sapply(features(fit), length)
#'
#' @author Martin Vincent
#' @importFrom sglOptim features
#' @export
features.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(features(object))
}

#' @title Nonzero parameters
#'
#' @description
#' Extracts the nonzero parameters for each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#'
#' @examples
#' data(SimData)
#'
#'
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl::fit(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero parameters of model 1, 10 and 25
#' parameters(fit)[c(1,10,25)]
#'
#' # count the number of nonzero parameters in each model
#' sapply(parameters(fit), sum)
#'
#' @author Martin Vincent
#' @importFrom sglOptim parameters
#' @export
parameters.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(parameters(object))
}


#' @title Extract feature statistics
#'
#' @description
#' Extracts the number of nonzero features (or group) in each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero features (or group) of the models.
#'
#' @author Martin Vincent
#' @importFrom sglOptim features_stat
#' @export
features_stat.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(features_stat(object, ...))
}


#' @title Extracting parameter statistics
#'
#' @description
#' Extracts the number of nonzero parameters in each model.
#'
#' @param object a msgl object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero parameters of the models.
#'
#' @author Martin Vincent
#' @importFrom sglOptim parameters_stat
#' @export
parameters_stat.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(parameters_stat(object, ...))
}


#' @title Number of models used for fitting
#' @description
#' Returns the number of models used for fitting.
#' Note that cv and subsampling objects does not containing any models even though nmod returns a positive number.
#'
#' @param object a msgl object
#' @param ... not used
#' @return the number of models in \code{object}
#'
#' @examples
#' data(SimData)
#'
#'
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl::fit(x, classes, alpha = .5, lambda = lambda)
#'
#' # the number of models
#' nmod(fit)
#'
#' @author Martin Vincent
#' @importFrom sglOptim nmod
#' @export
nmod.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(nmod(object, ...))
}

#' @title Index of best model
#'
#' @description
#' Returns the index of the best model, in terms of lowest error rate
#' @param object a msgl object
#' @param ... additional parameters (ignored)
#' @return index of the best model.
#'
#' @author Martin Vincent
#' @importFrom sglOptim best_model
#' @export
best_model.msgl <- function(object, ...) {
	class(object) <- "sgl" # Use std function
	return(best_model(object, "msgl", ...))
}


#' @title Extract the fitted models
#'
#' @description
#' Returns the fitted models, that is the estimated \eqn{\beta} matrices.
#'
#' @param object a msgl object
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of \eqn{\beta} matrices.
#'
#' @author Martin Vincent
#' @importFrom sglOptim models
#' @export
models.msgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(models(object, ...))
}

#' @title Nonzero coefficients
#' @description
#' This function returns the nonzero coefficients (that is the nonzero entries of the \eqn{beta} matrices)
#'
#' @param object a msgl object
#' @param index indices of the models
#' @param ... ignored
#' @return a list of length \code{length(index)} with nonzero coefficients of the models
#'
#' @examples
#' data(SimData)
#'
#'
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit <- msgl::fit(x, classes, alpha = .5, lambda = lambda)
#'
#' # the nonzero coefficients of the models 1, 10 and 20
#' coef(fit, index = c(1,10,20))
#'
#' @author Martin Vincent
#' @importFrom stats coef
#' @export
coef.msgl <- function(object, index = 1:nmod(object), ...) {
	class(object) <- "sgl" # Use std function
	return(coef(object, index = index, ...))
}


#' Print function for msgl
#'
#' This function will print some general information about the msgl object
#'
#' @param x msgl object
#' @param ... ignored
#'
#' @examples
#' data(SimData)
#'
#' ### Estimation
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 25, lambda.min = 0.075)
#' fit <- msgl::fit(x, classes, alpha = .5, lambda = lambda)
#'
#' # Print some information about the estimated models
#' fit
#'
#' ### Cross validation
#' fit.cv <- msgl::cv(x, classes, alpha = .5, lambda = lambda)
#'
#' # Print some information
#' fit.cv
#'
#' ### Subsampling
#' test <- list(1:20, 21:40)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#'
#' lambda <- msgl::lambda(x, classes, alpha = .5, d = 50, lambda.min = 0.05)
#' fit.sub <- msgl::subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test)
#'
#' # Print some information
#' fit.sub
#'
#' @author Martin Vincent
#' @importFrom sglOptim sgl_print
#' @method print msgl
#' @export
print.msgl <- function(x, ...) {
	sgl_print(x)
}
