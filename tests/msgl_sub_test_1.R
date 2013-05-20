library(msgl)

### Basic tests

data(SimData)
x <- sim.data$x
classes <- sim.data$classes

lambda <- msgl.lambda.seq(x, classes, alpha = .5, d = 100L, lambda.min = 0.01, standardize = TRUE)

test <- replicate(2, 1:20, simplify = FALSE)
train <- lapply(test, function(s) (1:length(classes))[-s])

fit.sub <- msgl.subsampling(x, classes, alpha = .5, lambda = lambda, training = train, test = test, max.threads = 2L)
if(!all(fit.sub$classes[[1]] == fit.sub$classes[[2]])) stop()
if(!all(colSums(fit.sub$classes[[2]] != classes[test[[2]]])[1:5*20] == c(12, 6, 5, 5, 4))) stop()


