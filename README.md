fastbackward
================

# Overview

**fastbackward** is a package that contains the `fastbackward()`
function. This function works similarly to backward elimination with the
`step()` function from the **stats** package; except, the
`fastbackward()` function makes use of a bounding algorithm to perform
backward elimination faster.

# How to install

It can be installed via the `install_github()` function from the
**devtools** package.

``` r
devtools::install_github("JacobSeedorff21/fastbackward")
```

# Usage

Here is a comparison of runtimes with the `step()` function from the
**stats** package and the `fastbackward()` function from the
**fastbackward** package. This comparison is based upon a randomly
generated logistic regression model with 1000 observations and 50
covariates.

``` r
# Loading in fastbackward
library(fastbackward)

# Defining function to generate datasets for logistic regression
LogisticSimul <- function(n, d, Bprob = .5, sd = 1, rho = 0.5){
  
  x <- MASS::mvrnorm(n, mu = rep(1, d), Sigma = diag(1 - rho, nrow = d, ncol = d) + 
                 matrix(rho, ncol = d, nrow = d))
  
  beta <- rnorm(d + 1, mean = 0, sd = sd) 
  
  beta[sample(2:length(beta), floor((length(beta) - 1) * Bprob))] = 0
  beta[beta != 0] <- beta[beta != 0] - mean(beta[beta != 0])
  
  p <- 1/(1 + exp(-x %*% beta[-1] - beta[1]))
  
  y <- rbinom(n, 1, p)
  
  df <- cbind(y, x) |> 
    as.data.frame()
  df
}

# Setting seed and creating dataset
set.seed(33391)
df <- LogisticSimul(1000, 50, .5, sd = 0.5)

# Fitting full logistic regression model
fullmodel <- glm(y ~ ., data = df, family = binomial(link = "logit"))

# Times
## Timing fast backward elimination
fastbackwardTime <- system.time(fastbackward1 <- fastbackward(fullmodel, trace = 0))
fastbackwardTime
```

    ##    user  system elapsed 
    ##    2.39    0.29    2.70

``` r
## Timing step function
stepTime <- system.time(BackwardStep <- step(fullmodel, direction = "backward", trace = 0))
stepTime
```

    ##    user  system elapsed 
    ##    8.41    0.97    9.42

For this logistic regression model, the fast backward elimination
algorithm from the **fastbackward** package was about 3.49 times faster
than step. The amount of speedup attained from the fast backward
elimination algorithm depends on the strength of association between the
covariates and the response variable. So, speedup will vary depending on
the specific problem.

### Checking results

``` r
# Checking if both methods give same results
all.equal(BackwardStep, fastbackward1)
```

    ## [1] TRUE

Hence, the two methods give the same results and the fast backward
elimination algorithm is faster than step.
