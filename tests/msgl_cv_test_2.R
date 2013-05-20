library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, sparse.data = TRUE, standardize = FALSE, max.threads = 2L, seed = 331L)
if(!all(colSums(fit.cv$classes != classes)[1:5*5] == c(50, 37, 29, 28, 25))) stop()
