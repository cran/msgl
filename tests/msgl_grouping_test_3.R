library(msgl)

### Tests grouping

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

### Define grouping

set.seed(100L)
grouping <- sample(1:100, replace = TRUE, size = 400)

## Lambda sequence
lambda <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.05, standardize = FALSE)

## Sparse group lasso

# Dense x
fit1a <- msgl(x, classes, grouping = grouping, alpha = .5, lambda = lambda, standardize = FALSE)
# (Forced) Sparse x
fit1b <- msgl(x, classes, grouping = grouping, alpha = .5, lambda = lambda, sparse.data = TRUE, standardize = FALSE)

if(max(abs(fit1a$beta[[100]]-fit1b$beta[[100]])) > 1e-10) stop()