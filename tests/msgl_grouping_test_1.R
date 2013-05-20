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

lambda1 <- msgl.lambda.seq(x, classes, grouping = grouping, alpha = .5, d = 100L, lambda.min = 0.05, sparse.data = TRUE, standardize = FALSE)
if(max(abs(lambda-lambda1)) != 0) stop()
