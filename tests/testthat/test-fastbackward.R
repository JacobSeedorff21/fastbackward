# Testing fast backward with lm
test_that("Testing fastbackward with lm", {
  library(fastbackward)
  set.seed(8621)
  x <- sapply(rep(0, 10), rnorm, n = 1000, simplify = TRUE)
  x <- cbind(1, x)
  beta <- rnorm(11)
  beta[3:9] <- 0
  y <- rnorm(n = 1000, mean = x %*% beta, sd = 2)
  Data <- cbind(y, x[,-1]) |>
    as.data.frame()

  ### Fitting full model
  fullmodel <- lm(y ~ ., data = Data)

  ### varying scale
  #### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  #### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  #### Checking results
  expect_equal(backward, fastbackward1)

  ### fixed scale
  #### step
  backward <- step(fullmodel, direction = "backward", trace = 0, scale = 2)

  #### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0, scale = 2)

  #### Checking results
  expect_equal(backward, fastbackward1)

  ### Fitting full model with interactions
  Data <- Data[, 1:6]
  fullmodel <- lm(y ~ .^2, data = Data)

  ### varying scale
  #### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  #### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  #### Checking results
  expect_equal(backward, fastbackward1)

  ### fixed scale
  #### step
  backward <- step(fullmodel, direction = "backward", trace = 0, scale = 2)

  #### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0, scale = 2)

  #### Checking results
  expect_equal(backward, fastbackward1)
})

## Testing different ks
test_that("Testing fastbackward with lm", {
  library(fastbackward)
  set.seed(8621)
  x <- sapply(rep(0, 10), rnorm, n = 1000, simplify = TRUE)
  x <- cbind(1, x)
  beta <- rnorm(11)
  beta[3:9] <- 0
  y <- rnorm(n = 1000, mean = x %*% beta, sd = 2)
  Data <- cbind(y, x[,-1]) |>
    as.data.frame()

  ### Fitting full model
  fullmodel <- lm(y ~ ., data = Data)

  ### varying scale
  #### step
  backward0 <- step(fullmodel, direction = "backward", trace = 0, k = 0)
  backward10 <- step(fullmodel, direction = "backward", trace = 0, k = 10)
  backwardlog <- step(fullmodel, direction = "backward", trace = 0, k = log(1000))

  #### fast backward
  fastbackward0 <- fastbackward(fullmodel, trace = 0, k = 0)
  fastbackward10 <- fastbackward(fullmodel, trace = 0, k = 10)
  fastbackwardlog <- fastbackward(fullmodel, trace = 0, k = log(1000))

  #### Checking results
  expect_equal(backward0, fastbackward0)
  expect_equal(backward10, fastbackward10)
  expect_equal(backwardlog, fastbackwardlog)

  ### fixed scale
  #### step
  backward0 <- step(fullmodel, direction = "backward", trace = 0, k = 0)
  backward10 <- step(fullmodel, direction = "backward", trace = 0, k = 10)
  backwardlog <- step(fullmodel, direction = "backward", trace = 0, k = log(1000))

  #### fast backward
  fastbackward0 <- fastbackward(fullmodel, trace = 0, k = 0)
  fastbackward10 <- fastbackward(fullmodel, trace = 0, k = 10)
  fastbackwardlog <- fastbackward(fullmodel, trace = 0, k = log(1000))

  #### Checking results
  expect_equal(backward0, fastbackward0)
  expect_equal(backward10, fastbackward10)
  expect_equal(backwardlog, fastbackwardlog)

  ### Fitting full model with interactions
  Data <- Data[, 1:6]
  fullmodel <- lm(y ~ .^2, data = Data)

  ### varying scale
  #### step
  backward0 <- step(fullmodel, direction = "backward", trace = 0, k = 0)
  backward10 <- step(fullmodel, direction = "backward", trace = 0, k = 10)
  backwardlog <- step(fullmodel, direction = "backward", trace = 0, k = log(1000))

  #### fast backward
  fastbackward0 <- fastbackward(fullmodel, trace = 0, k = 0)
  fastbackward10 <- fastbackward(fullmodel, trace = 0, k = 10)
  fastbackwardlog <- fastbackward(fullmodel, trace = 0, k = log(1000))

  #### Checking results
  expect_equal(backward0, fastbackward0)
  expect_equal(backward10, fastbackward10)
  expect_equal(backwardlog, fastbackwardlog)

  ### fixed scale
  #### step
  backward0 <- step(fullmodel, direction = "backward", trace = 0, k = 0)
  backward10 <- step(fullmodel, direction = "backward", trace = 0, k = 10)
  backwardlog <- step(fullmodel, direction = "backward", trace = 0, k = log(1000))

  #### fast backward
  fastbackward0 <- fastbackward(fullmodel, trace = 0, k = 0)
  fastbackward10 <- fastbackward(fullmodel, trace = 0, k = 10)
  fastbackwardlog <- fastbackward(fullmodel, trace = 0, k = log(1000))

  #### Checking results
  expect_equal(backward0, fastbackward0)
  expect_equal(backward10, fastbackward10)
  expect_equal(backwardlog, fastbackwardlog)
})

# Testing fast backward with glm
## Gamma regression
test_that("Testing fast backward with gamma regression", {
  library(fastbackward)
  set.seed(8623)

  ### Making sure x and beta are positive to use inverse link
  x <- sapply(rep(1, 10), rexp, n = 1000, simplify = TRUE)
  x <- cbind(1, x)
  beta <- rexp(11, rate = 15)
  beta[3:6] <- 0
  y <- rgamma(n = 1000, shape = 1, scale = 1 / (x %*% beta))
  Data <- cbind(y, x[,-1]) |>
    as.data.frame()

  ### Fitting full model
  fullmodel <- glm(y ~ ., data = Data, family = Gamma(link = "inverse"))

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)

  ### Fitting full model with interactions
  Data <- Data[, 1:6]
  fullmodel <- glm(y ~ .^2, data = Data, family = Gamma(link = "inverse"))

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)
})

## Logistic regression
test_that("Testing fast backward with logistic regression", {
  library(fastbackward)
  set.seed(8621)
  x <- sapply(rep(0, 10), rnorm, n = 1000, simplify = TRUE)
  x <- cbind(1, x)
  beta <- rnorm(11)
  y <- rbinom(n = 1000, size = 1, p = plogis(x %*% beta + rnorm(1000, sd = 3)))
  Data <- cbind(y, x[,-1]) |>
    as.data.frame()

  ### Fitting full model
  fullmodel <- glm(y ~ ., data = Data, family = binomial)

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)

  ### Fitting full model with interactions
  Data <- Data[, 1:6]
  fullmodel <- glm(y ~ .^2, data = Data, family = binomial)

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)
})

## Poisson
test_that("Testing fast backward with poisson regression", {
  library(fastbackward)
  set.seed(8621)
  x <- sapply(rep(0, 10), rnorm, n = 1000, simplify = TRUE)
  x <- cbind(1, x)
  beta <- rnorm(11, sd = 0.1)
  beta[c(2, 5, 6, 7, 10)] <- 0
  y <- rpois(n = 1000, lambda = exp(x %*% beta))
  Data <- cbind(y, x[,-1]) |>
    as.data.frame()

  ### Fitting full model
  fullmodel <- glm(y ~ ., data = Data, family = poisson)

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)

  ### Fitting full model with interactions
  Data <- Data[, 1:6]
  fullmodel <- glm(y ~ .^2, data = Data, family = poisson)

  ### step
  backward <- step(fullmodel, direction = "backward", trace = 0)

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0)

  ### Checking results
  expect_equal(backward, fastbackward1)
})
