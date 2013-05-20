library(msgl)

### Test r out stream

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

config.verbose <- sgl.algorithm.config(verbose = TRUE)

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

# msgl.cv
fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, standardize = FALSE, max.threads = 2L, seed = 331L, algorithm.config = config.verbose)
