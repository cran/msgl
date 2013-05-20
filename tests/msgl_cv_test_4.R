library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = FALSE)

fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, standardize = FALSE, max.threads = 2L, seed = 331L)
if(!all(colSums(fit.cv$classes != classes)[1:5*5] == c(67, 48, 33, 28, 25))) stop()