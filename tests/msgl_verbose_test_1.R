library(msgl)

### Test r out stream

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

config.verbose <- sgl.algorithm.config(verbose = TRUE)

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

# msgl
fit <- msgl(x, classes, alpha = 0, lambda = lambda, standardize = TRUE, algorithm.config = config.verbose)
