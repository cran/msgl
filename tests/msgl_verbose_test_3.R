library(msgl)

### Test r out stream

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

config.verbose <- sgl.algorithm.config(verbose = TRUE)

## Lambda sequence

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 25L, lambda.min = 0.05, standardize = TRUE)

# msgl.subsampling
test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test, max.threads = 2L, algorithm.config = config.verbose)
