# Testing scope
test_that("Testing scope argument", {
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
  backward1 <- step(fullmodel, direction = "backward", trace = 0, scope = formula("y ~ 1"))
  backward2 <- step(fullmodel, direction = "backward", trace = 0, scope = list("lower" = formula("y ~ V3 + V4 + V6")))
  backward3 <- step(fullmodel, direction = "backward", trace = 0, scope = list("lower" = formula("y ~ .")))
  backward4 <- step(fullmodel, direction = "backward", trace = 0, scope = list("lower" = formula("y ~ . - V2 - V4 - V11")))

  ### fast backward
  fastbackward1 <- fastbackward(fullmodel, trace = 0, scope = formula("y ~ 1"))
  fastbackward2 <- fastbackward(fullmodel, trace = 0, scope = formula("y ~ V3 + V4 + V6"))
  fastbackward3 <- fastbackward(fullmodel, trace = 0, scope = formula("y ~ ."))
  fastbackward4 <- fastbackward(fullmodel, trace = 0, scope = formula("y ~ . - V2 - V4 - V11"))

  ### Checking results
  expect_equal(backward1, fastbackward1)
  expect_equal(backward2, fastbackward2)
  expect_equal(backward3, fastbackward3)
  expect_equal(backward4, fastbackward4)
})

# Testing bad inputs to scope
test_that("Testing bad inputs to scope argument", {
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

  ### fast backward
  expect_error(fastbackward1 <- fastbackward(fullmodel, trace = 0, scope = "apple"))
  expect_error(fastbackward2 <- fastbackward(fullmodel, trace = 0, scope = list("lower" =  formula("y ~ V3 + V4 + V6"))))
  expect_error(fastbackward3 <- fastbackward(fullmodel, trace = 0, scope = 6))
  expect_error(fastbackward4 <- fastbackward(fullmodel, trace = 0, scope = NA))
})
