library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

set.seed(100L)

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

fit.cv <- msgl.cv(x, classes, alpha = .5, fold = 11L, lambda = lambda, standardize = TRUE, max.threads = 2L, seed = 331L)
if(!all(colSums(fit.cv$classes != classes)[1:5*5] == c(98, 70, 53, 49, 42))) stop()
