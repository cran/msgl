library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = TRUE)

fit.cv <- msgl.cv(x, classes, alpha = .5, lambda = lambda, standardize = TRUE, max.threads = 2L, seed = 331L)
if(!all(colSums(fit.cv$classes != classes)[1:5*20] == c(36, 25, 18, 15, 16))) stop()