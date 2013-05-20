library(msgl)

### Test r out stream

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

config.verbose <- sgl.algorithm.config(verbose = TRUE)

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = TRUE)

# msgl
fit <- msgl(x, classes, alpha = 0, lambda = lambda, standardize = TRUE, algorithm.config = config.verbose)
