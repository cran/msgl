# TODO: Add comment
# 
# Author: martin
###############################################################################

#' Fit a multinomial sparse group lasso regularization path. 
#'
#' For a classification problem with  \eqn{K} classes and \eqn{p} covariates dived into \eqn{m} groups.
#' A sequence of minimizers (one for each lambda given in the \code{lambda} argument) of 
#' \deqn{\hat R(\beta) + \lambda \left( (1-\alpha) \sum_{J=1}^m \gamma_J \|\beta^{(J)}\|_2 + \alpha \sum_{i=1}^{n} \xi_i |\beta_i| \right)}
#' where \eqn{\hat R} is the weighted empirical log-likelihood risk of the multinomial regression model.
#' The vector \eqn{\beta^{(J)}} denotes the parameters associated with the \eqn{J}'th group of covariates
#' (default is one covariate per group, hence the default dimension of \eqn{\beta^{(J)}} is \eqn{K}). 
#' The group weights \eqn{\gamma \in [0,\infty)^m} and the parameter weights \eqn{\xi = (\xi^{(1)},\dots, \xi^{(m)}) \in [0,\infty)^n} 
#' with \eqn{\xi^{(1)}\in [0,\infty)^{n_1},\dots, \xi^{(m)} \in [0,\infty)^{n_m}}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param return the indices of lambda values for which to return a the fitted parameters.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{beta}{the fitted parameters -- a list of length \code{length(lambda)} with each entry a matrix of size \eqn{K\times (p+1)} holding the fitted parameters}
#' \item{loss}{the values of the loss function}
#' \item{objective}{the values of the objective function (i.e. loss + penalty)}
#' \item{lambda}{the lambda values used}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' fit <- msgl(x, classes, lambda = lambda)
#' fit$beta[[10]] #model with lambda = lambda[10] 
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_wb_basic r_msgl_wb_sp_basic
#' @import Matrix
msgl <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, return = 1:length(lambda), sparse.data = FALSE, algorithm.config = sgl.standard.config) {
	
	classes <- factor(classes)
	nclasses <- length(levels(classes))
	
	if(!is.null(grouping)) {
		
		grouping <- factor(grouping)
		group.order <- order(grouping)
		
		#Reorder 
		x <- x[,group.order]
		
		if(is.null(groupWeights)) {
			groupWeights <- c(0, sqrt(nclasses*table(grouping))) #FIXME intercept weight
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights <- matrix(1, nrow = nclasses, ncol = ncol(x) + 1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- nclasses*c(1L, as.integer(table(grouping)))
		
	} else {	
		
		if(is.null(groupWeights)) {
			groupWeights = c(0,rep(sqrt(nclasses), ncol(x)))
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights = matrix(1, nrow = nclasses, ncol = ncol(x)+1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- rep(nclasses, ncol(x)+1)
	}
	
	#TODO check that return is valid
	return <- as.integer(sort(unique(return))) - 1L
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	classes.numeric <- as.integer(factor(classes))-1L
	
	if(standardize) {
		
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- if(sparse.data) rep(0, length(x.scale)) else attr(x, "scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	# Run msgl algorithm
	
	if(sparse.data) {
		
		m <- as(x, "CsparseMatrix")
		m <- list(dim(m), m@p, m@i, m@x)
		
		res <- .Call(r_msgl_wb_sp_basic, m, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, return, algorithm.config)
	} else {
		res <- .Call(r_msgl_wb_basic, x, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, return, algorithm.config)
	}
	
	# Dim names
	feature.names <- if(!is.null(colnames(x))) colnames(x) else 1:dim(x)[2]
	
	if(is.list(classes)) {
		
		class.names <- unlist(lapply(classes, function(x) levels(factor(x))))
		
	} else {
		
		class.names <- levels(factor(classes))
		
	}
	
	# Create R sparse matrix
	res$beta <- lapply(1:length(res$beta), function(i) sparseMatrix(p = res$beta[[i]][[2]], i = res$beta[[i]][[3]], x = res$beta[[i]][[4]], dims = res$beta[[i]][[1]], dimnames = list(class.names, c("Intercept", feature.names)), index1 = FALSE))
	
	# Convert beta back to the org scale
	if(standardize) {
		
		res$beta <- .to_org_scale(beta = res$beta, x.scale = x.scale, x.center = x.center)
		
	}
	
	if(!is.null(grouping)) {
		
		# Restore org order
		res$beta <- lapply(res$beta, function(beta.matrix) beta.matrix[, c(1,1+order(group.order))])
	}
	
	#TODO name
	res$classes.prior <- classes
	
	class(res) <- "msgl"
	return(res)
}

#' Computes a lambda sequence for the regularization path
#' 
#' Computes a decreasing lambda sequence of length \code{d}.
#' The sequence ranges from a data determined maximal lambda \eqn{\lambda_\textrm{max}} to the user inputed \code{lambda.min}.
#'
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param d the length of lambda sequence
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda.min the smallest lambda value in the computed sequence. 
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param algorithm.config the algorithm configuration to be used. 
#' @return a vector of length \code{d} containing the compute lambda sequence.
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_wb_sp_lambda_seq r_msgl_wb_lambda_seq
msgl.lambda.seq <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, d = 100L, standardize = TRUE, lambda.min, sparse.data = FALSE, algorithm.config = sgl.standard.config) {
	
	classes <- factor(classes)
	nclasses <- length(levels(classes))
	
	if(!is.null(grouping)) {
		
		grouping <- factor(grouping)
		group.order <- order(grouping)
		
		#Reorder 
		x <- x[,group.order]
		
		if(is.null(groupWeights)) {
			groupWeights <- c(0, sqrt(nclasses*table(grouping))) #FIXME intercept weight
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights <- matrix(1, nrow = nclasses, ncol = ncol(x) + 1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- nclasses*c(1L, as.integer(table(grouping)))
		
	} else {	
		
		if(is.null(groupWeights)) {
			groupWeights = c(0,rep(sqrt(nclasses), ncol(x)))
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights = matrix(1, nrow = nclasses, ncol = ncol(x)+1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- rep(nclasses, ncol(x)+1)
	}
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	classes.numeric <- as.integer(factor(classes))-1L
	
	if(standardize) {
		
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- attr(x, "scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	if(sparse.data) {
		
		m <- as(x, "CsparseMatrix")
		m <- list(dim(m), m@p, m@i, m@x)
		
		res <- .Call(r_msgl_wb_sp_lambda_seq, m, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, d, lambda.min, algorithm.config)
	} else {
		res <- .Call(r_msgl_wb_lambda_seq, x, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, d, lambda.min, algorithm.config)
	}
	
	return(res)
}

#' Multinomial sparse group lasso cross validation using multiple possessors 
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param fold the fold of the cross validation, an integer larger than \eqn{1} and less than \eqn{N+1}. Ignored if \code{cv.indices != NULL}.
#' If \code{fold}\eqn{\le}\code{max(table(classes))} then the data will be split into \code{fold} disjoint subsets keeping the ration of classes approximately equal.
#' Otherwise the data will be split into \code{fold} disjoint subsets without keeping the ration fixed.
#' @param cv.indices a list of indices of a cross validation splitting. 
#' If \code{cv.indices = NULL} then a random splitting will be generated using the \code{fold} argument.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param max.threads the maximal number of threads to be used
#' @param seed the seed used for generating the random cross validation splitting, only used if \code{fold}\eqn{\le}\code{max(table(classes))}. 
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the linear predictors.}
#' \item{response}{the estimated probabilities - a list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the probabilities.}
#' \item{classes}{the estimated classes - a matrix of size \eqn{N \times d} with \eqn{d=}\code{length(lambda)}.}
#' \item{cv.indices}{the cross validation splitting used.}
#' \item{features}{average number of features used in the models.}
#' \item{parameters}{average number of parameters used in the models.}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 50L, lambda.min = 0.03)
#' fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda)
#' 
#' # Missclassification count
#' colSums(fit.cv$classes != classes)
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_wb_cv r_msgl_wb_sp_cv
msgl.cv <- function(x, classes, sampleWeights = NULL, grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, fold = 10L, cv.indices = list(), sparse.data = FALSE, max.threads = 2L, seed = 331L, algorithm.config = sgl.standard.config) {
	
	classes <- factor(classes)
	nclasses <- length(levels(classes))
	
	if(is.null(sampleWeights)) {
		if(length(cv.indices) == 0) {
			sampleWeights <- rep(fold/(length(classes)*(fold-1)), length(classes))
		} else {
			n_train <- sapply(cv.indices, function(x) length(classes)-length(x))
			sampleWeights <- rep(1/mean(n_train), length(classes))
		}
	}
	
	if(!is.null(grouping)) {
		
		grouping <- factor(grouping)
		group.order <- order(grouping)
		
		#Reorder 
		x <- x[,group.order]
		
		if(is.null(groupWeights)) {
			groupWeights <- c(0, sqrt(nclasses*table(grouping))) #FIXME intercept weight
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights <- matrix(1, nrow = nclasses, ncol = ncol(x) + 1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- nclasses*c(1L, as.integer(table(grouping)))
		
	} else {	
		
		if(is.null(groupWeights)) {
			groupWeights = c(0,rep(sqrt(nclasses), ncol(x)))
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights = matrix(1, nrow = nclasses, ncol = ncol(x)+1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- rep(nclasses, ncol(x)+1)
	}
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	classes.numeric <- as.integer(factor(classes))-1L
	
	if(standardize) {
		
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- attr(x, "scaled:center")
	}
	
	#TODO domain check
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	
	#Dimnames
	class.names <- levels(factor(classes))
	dim.names <-  list(class.names, rownames(x))
	
	if(length(cv.indices) == 0) {
		
		use.cv.indices <- FALSE
		
		# Check fold
		if(fold < 2) {
			stop("fold must be equal to or larger than 2")
		}
		
		if(fold > length(classes)) {
			stop("fold must be equal to or less than the number of samples")
		}
		
		if(fold > max(table(classes))) {
			# use random sample indices
			use.cv.indices <- TRUE
			cv.indices <- split(sample(0:(length(classes))-1L), 1:fold)
		}
		
	} else {
		
		cv.indices <- lapply(cv.indices, function(x) as.integer(x-1))
		use.cv.indices <- TRUE
	}
	
	if(sparse.data) {
		
		m <- as(x, "CsparseMatrix")
		m <- list(dim(m), m@p, m@i, m@x)	
		
		res <- .Call(r_msgl_wb_sp_cv, m, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, use.cv.indices, max.threads, seed, algorithm.config)				
	} else {
		res <- .Call(r_msgl_wb_cv, x, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, fold, cv.indices, use.cv.indices, max.threads, seed, algorithm.config)				
	}
	
	#Set class names
	rownames(res$classes) <- dim.names[[2]]
	
	if(!is.null(dim.names[[1]])) {
		res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
	}
	
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <-dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <-dim.names; m})
	
	
	class(res) <- "msgl"
	return(res)
}

#' Multinomial sparse group lasso generic subsampling procedure
#'
#' Support the use of multiple processors.
#' 
#' @param x design matrix, matrix of size \eqn{N \times p}.
#' @param classes classes, factor of length \eqn{N}.
#' @param sampleWeights sample weights, a vector of length \eqn{N}.
#' @param grouping grouping of covariates, a vector of length \eqn{p}. Each element of the vector specifying the group of the covariate. 
#' @param groupWeights the group weights, a vector of length \eqn{m+1} (the number of groups). 
#' The first element of the vector is the intercept weight. 
#' If \code{groupWeights = NULL} default weights will be used.
#' Default weights are 0 for the intercept and \deqn{\sqrt{K\cdot\textrm{number of covariates in the group}}} for all other weights.
#' @param parameterWeights a matrix of size \eqn{K \times (p+1)}. 
#' The first column of the matrix is the intercept weights.
#' Default weights are is 0 for the intercept weights and 1 for all other weights.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param standardize if TRUE the covariates are standardize before fitting the model. The model parameters are returned in the original scale. 
#' @param lambda the lambda sequence for the regularization path.
#' @param training a list of training samples, each item of the list corresponding to a subsample.
#' Each item in the list must be a vector with the indices of the training samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{test} list.  
#' @param test a list of test samples, each item of the list corresponding to a subsample.
#' Each item in the list must be vector with the indices of the test samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{training} list.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param max.threads the maximal number of threads to be used
#' @param algorithm.config the algorithm configuration to be used. 
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(test)} with each element of the list another list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the linear predictors.}
#' \item{response}{the estimated probabilities -- a list of length \code{length(test)} with each element of the list another list of length \code{length(lambda)} one item for each lambda value, with each item a matrix of size \eqn{K \times N} containing the probabilities.}
#' \item{classes}{the estimated classes -- a list of length \code{length(test)} with each element of the list a matrix of size \eqn{N \times d} with \eqn{d=}\code{length(lambda)}.}
#' \item{features}{number of features used in the models.}
#' \item{parameters}{number of parameters used in the models.}
#' @examples 
#' data(SimData)
#' x <- sim.data$x
#' classes <- sim.data$classes
#' lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.03)
#'
#' test <- replicate(5, sample(1:length(classes))[1:20], simplify = FALSE)
#' train <- lapply(test, function(s) (1:length(classes))[-s])
#' 
#' fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, 
#'  training = train, test = test)
#' 
#' # Missclassification count of second subsample
#' colSums(fit.sub$classes[[2]] != classes[test[[2]]])
#' @author Martin Vincent
#' @export
#' @useDynLib msgl r_msgl_wb_subsampling r_msgl_wb_sp_subsampling
msgl.subsampling <- function(x, classes, sampleWeights = rep(1/length(classes), length(classes)), grouping = NULL, groupWeights = NULL, parameterWeights = NULL, alpha = 0.5, standardize = TRUE, 
		lambda, training, test, sparse.data = FALSE, max.threads = 2L, algorithm.config = sgl.standard.config) {
	
	classes <- factor(classes)
	nclasses <- length(levels(classes))
	
	if(!is.null(grouping)) {
		
		grouping <- factor(grouping)
		group.order <- order(grouping)
		
		#Reorder 
		x <- x[,group.order]
		
		if(is.null(groupWeights)) {
			groupWeights <- c(0, sqrt(nclasses*table(grouping))) #FIXME intercept weight
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights <- matrix(1, nrow = nclasses, ncol = ncol(x) + 1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- nclasses*c(1L, as.integer(table(grouping)))
		
	} else {	
		
		if(is.null(groupWeights)) {
			groupWeights = c(0,rep(sqrt(nclasses), ncol(x)))
		}
		
		if(is.null(parameterWeights)) {
			parameterWeights = matrix(1, nrow = nclasses, ncol = ncol(x)+1)
			parameterWeights[,1] <- 0
		}
		
		#Compute block dim
		block.dim <- rep(nclasses, ncol(x)+1)
	}
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	classes.numeric <- as.integer(factor(classes))-1L
	
	if(standardize) {
		
		x <- scale(x, if(sparse.data) FALSE else TRUE, TRUE)
		x.scale <- attr(x, "scaled:scale")
		x.center <- attr(x, "scaled:center")
	}
	
	#TODO domain check
	
	if(sum(is.na(x)) > 0) {
		warning("Replacing NA values in x with 0")
		x[is.na(x)] <- 0
	}
	
	#Dimnames
	
	class.names <- levels(factor(classes))
	dim.names <-  list(class.names, rownames(x))
	
	trainingSamples.0 <- lapply(training, function(x) as.integer(x - 1))
	testSamples.0 <- lapply(test, function(x) as.integer(x - 1))
	
	if(sparse.data) {
		
		m <- as(x, "CsparseMatrix")
		m <- list(dim(m), m@p, m@i, m@x)
		
		res <- .Call(r_msgl_wb_sp_subsampling, m, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, trainingSamples.0, testSamples.0, max.threads, algorithm.config)	
	} else {
		res <- .Call(r_msgl_wb_subsampling, x, classes.numeric, sampleWeights, block.dim, groupWeights, parameterWeights, alpha, lambda, trainingSamples.0, testSamples.0, max.threads, algorithm.config)	
	}

	for(i in 1:length(test)) {
		res$classes[[i]] <- res$classes[[i]]+1
	}
		
	for(i in 1:length(test)) {
		
		#Set class names
		rownames(res$classes[[i]]) <- dim.names[[2]][test[[i]]]
		
		if(!is.null(dim.names[[1]])) {
			res$classes[[i]] <- apply(X = res$classes[[i]], MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x])
		}

		res$link[[i]] <- lapply(X = res$link[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][test[[i]]]); m})
		res$response[[i]] <- lapply(X = res$response[[i]], FUN = function(m) {dimnames(m) <- list(dim.names[[1]], dim.names[[2]][test[[i]]]); m})
	}
	
	class(res) <- "msgl"
	return(res)
}


#' Create a new algorithm configuration
#'
#' With the exception of \code{verbose} it is not recommended to change any of the default values.
#' 
#' @param tolerance_penalized_main_equation_loop tolerance threshold.
#' @param tolerance_penalized_inner_loop_alpha tolerance threshold.
#' @param tolerance_penalized_inner_loop_beta tolerance threshold.
#' @param tolerance_penalized_middel_loop_alpha tolerance threshold.
#' @param tolerance_penalized_outer_loop_alpha tolerance threshold.
#' @param tolerance_penalized_outer_loop_beta tolerance threshold.
#' @param tolerance_penalized_outer_loop_gamma tolerance threshold.
#' @param use_bound_optimization if \code{TRUE} hessian bound check will be used.
#' @param use_stepsize_optimization_in_penalizeed_loop if \code{TRUE} step-size optimization will be used.
#' @param stepsize_opt_penalized_initial_t initial step-size.
#' @param stepsize_opt_penalized_a step-size optimization parameter.
#' @param stepsize_opt_penalized_b step-size optimization parameter.
#' @param verbose If \code{TRUE} some information, regarding the status of the algorithm, will be printed in the R terminal.
#' @return A configuration.
#' @examples
#' config.verbose <- sgl.algorithm.config(verbose = TRUE)
#' @author Martin Vincent
#' @export
sgl.algorithm.config <- function(tolerance_penalized_main_equation_loop = 1e-10, 
		tolerance_penalized_inner_loop_alpha = 1e-4, 
		tolerance_penalized_inner_loop_beta = 1, 
		tolerance_penalized_middel_loop_alpha = 0.01, 
		tolerance_penalized_outer_loop_alpha = 0.01, 
		tolerance_penalized_outer_loop_beta = 0, 
		tolerance_penalized_outer_loop_gamma = 1e-5, 
		use_bound_optimization = TRUE, 
		use_stepsize_optimization_in_penalizeed_loop = TRUE, 
		stepsize_opt_penalized_initial_t = 1,
		stepsize_opt_penalized_a = 0.1, 
		stepsize_opt_penalized_b = 0.1, 
		verbose = FALSE) {
	
	config <- list()
	
	config$tolerance_penalized_main_equation_loop <- tolerance_penalized_main_equation_loop
	
	config$tolerance_penalized_inner_loop_alpha <- tolerance_penalized_inner_loop_alpha
	config$tolerance_penalized_inner_loop_beta <- tolerance_penalized_inner_loop_beta
	
	config$tolerance_penalized_middel_loop_alpha <- tolerance_penalized_middel_loop_alpha
	
	config$tolerance_penalized_outer_loop_alpha <- tolerance_penalized_outer_loop_alpha
	config$tolerance_penalized_outer_loop_beta <- tolerance_penalized_outer_loop_beta	
	config$tolerance_penalized_outer_loop_gamma <- tolerance_penalized_outer_loop_gamma	
	
	config$use_bound_optimization <- use_bound_optimization
	
	config$use_stepsize_optimization_in_penalizeed_loop <- use_stepsize_optimization_in_penalizeed_loop
	config$stepsize_opt_penalized_initial_t <- stepsize_opt_penalized_initial_t
	config$stepsize_opt_penalized_a <- stepsize_opt_penalized_a
	config$stepsize_opt_penalized_b <- stepsize_opt_penalized_b
	
	config$verbose <- verbose
	
	return(config)
}

#' Standard algorithm configuration
#'
#' \code{sgl.standard.config <- sgl.algorithm.config()}
#' 
#' @author Martin Vicnet
#' @export
sgl.standard.config <- sgl.algorithm.config();

.to_org_scale <- function(beta, x.scale, x.center) {
	for(l in 1:length(beta)) {
		
		beta.org <- t(t(beta[[l]])*c(1,1/x.scale))
		beta.org[,1] <- beta.org[,1] - rowSums(t(t(beta[[l]][,-1])*(x.center/x.scale)))
		
		beta[[l]] <- beta.org
	}
	
	return(beta)
}

.to_std_scale <- function(beta, x.scale, x.center) {
	
	for(l in 1:length(beta)) {
		beta.std <- t(t(beta[[l]])*c(1, x.scale))
		beta.std[,1] <- beta.std[,1] + rowSums(t(t(beta[[l]][,-1])*(x.center)))
		
		beta[[l]] <- beta.std
	}
	
	return(beta)
}

#' Predict
#' 
#' Computes the linear predictors, the estimated probabilities and the estimated classes for a new data set.
#'
#' @param object an object of class msgl, produced with \code{msgl}.
#' @param x a data matrix of size \eqn{N_\textrm{new} \times p}.
#' @param sparse.data if TRUE \code{x} will be treated as sparse, if \code{x} is a sparse matrix it will be treated as sparse by default.
#' @param ... ignored.
#' @return 
#' \item{link}{the linear predictors -- a list of length \code{length(fit$beta)} one item for each model, with each item a matrix of size \eqn{K \times N_\textrm{new}} containing the linear predictors.}
#' \item{response}{the estimated probabilities -- a list of length \code{length(fit$beta)} one item for each model, with each item a matrix of size \eqn{K \times N_\textrm{new}} containing the probabilities.}
#' \item{classes}{the estimated classes -- a matrix of size \eqn{N_\textrm{new} \times d} with \eqn{d=}\code{length(fit$beta)}.}
#' @author Martin Vincent
#' @method predict msgl
#' @S3method predict msgl
#' @export
#' @useDynLib msgl r_msgl_predict r_msgl_sparse_predict
predict.msgl <- function(object, x, sparse.data = FALSE, ...) {
	
	if(!sparse.data) {
		sparse.data <- is(x, "sparseMatrix")
	}
	
	res <- list()
	
	if("beta" %in% names(object)) {
		res <- .predict.msgl(object$beta, x, sparse.data)	
	} else  {
		stop("No models found - missing beta")
	}
	
	class(res) <- "msgl"
	return(res)
}

.predict.msgl <- function(beta, x, sparse.data = FALSE) {
	#Save dimnames
	dim.names <-  list(rownames(beta[[1]]), rownames(x))
	
	beta <- lapply(X = beta, FUN = function(m) as(m, "CsparseMatrix"))
	beta <- lapply(X = beta, FUN = function(m) list(dim(m), m@p, m@i, m@x))
	
	if(sparse.data) {
		
		m <- as(x, "CsparseMatrix")
		m <- list(dim(m), m@p, m@i, m@x)
		
		res <- .Call(r_msgl_sparse_predict, m, beta)
		
	} else {
		
		res <- .Call(r_msgl_predict, x, beta)
	}
	
	#Set class names
	rownames(res$classes) <- dim.names[[2]]
	
	if(!is.null(dim.names[[1]])) {
		res$classes <- apply(X = res$classes, MARGIN = c(1,2), FUN = function(x) dim.names[[1]][x+1])
	}
	
	#Set dimnames ect
	
	res$link <- lapply(X = res$link, FUN = function(m) {dimnames(m) <-dim.names; m})
	res$response <- lapply(X = res$response, FUN = function(m) {dimnames(m) <-dim.names; m})
	
	
	class(res) <- "msgl"
	return(res)
}

#' Simulated data set
#'
#' The use of this data set is only intended for testing and examples.
#' The data set contains 100 simulated samples grouped into 10 classes.
#' For each sample 400 covariates have been simulated.
#'
#' @name sim.data
#' @docType data
#' @keywords data
NULL

