# ---- Box-Cox Tests ----

test_that("Box-Cox transformation and inverse are consistent", {
  y <- c(1, 2, 5, 10)
  lambda <- 0.5
  f <- tstransform(method = "box-cox", lambda = lambda)
  transformed <- f$transform(y)
  recovered <- f$inverse(transformed, lambda = lambda)
  expect_equal(as.numeric(recovered), y, tolerance = 1e-6)
})

test_that("Box-Cox handles lambda = 0 correctly", {
  y <- c(1, 2, 5, 10)
  f <- tstransform(method = "box-cox", lambda = 0)
  transformed <- f$transform(y)
  recovered <- f$inverse(transformed, lambda = 0)
  expect_equal(as.numeric(recovered), y, tolerance = 1e-6)
})

test_that("Box-Cox errors on negative values for lambda < 0", {
  y <- c(-1, 0, 2, 5)
  f <- tstransform(method = "box-cox", lambda = -0.5)
  expect_error(f$transform(y))
})


# ---- Logit Tests ----
test_that("Logit transformation and inverse are consistent", {
  x <- seq(0.1, 0.9, length.out = 10)
  f <- tstransform(method = "logit", lower = 0, upper = 1)
  transformed <- f$transform(x)
  recovered <- f$inverse(transformed)
  expect_equal(as.numeric(recovered), x, tolerance = 1e-6)
})

test_that("Logit errors on values outside bounds", {
  f <- tstransform(method = "logit", lower = 0, upper = 1)
  expect_error(f$transform(c(-0.1, 1.1)))
})

# ---- Softplus Logit Tests ----
test_that("Softplus Logit transformation and inverse are consistent", {
  x <- seq(0.1, 0.9, length.out = 10)
  f <- tstransform(method = "softplus-logit", lower = 0, upper = 1)
  transformed <- f$transform(x)
  recovered <- f$inverse(transformed)
  expect_equal(as.numeric(recovered), x, tolerance = 1e-6)
})

test_that("Softplus Logit errors on values outside bounds", {
  f <- tstransform(method = "softplus-logit", lower = 0, upper = 1)
  expect_error(f$transform(c(-0.1, 1.1)))
})


# ---- Sigmoid Tests ----
test_that("Sigmoid transformation and inverse are consistent", {
  x <- seq(-5, 5, length.out = 10)
  f <- tstransform(method = "sigmoid", lower = 0, upper = 1)
  transformed <- f$transform(x)
  recovered <- f$inverse(transformed)
  expect_equal(as.numeric(recovered), x, tolerance = 1e-6)
})

test_that("Sigmoid produces values strictly between lower and upper", {
  x <- seq(-5, 5, length.out = 10)
  f <- tstransform(method = "sigmoid", lower = 0, upper = 1)
  transformed <- f$transform(x)
  expect_true(all(transformed > 0 & transformed < 1))
})

test_that("Sigmoid inverse errors on values outside bounds", {
  f <- tstransform(method = "sigmoid", lower = 0, upper = 1)
  expect_error(f$inverse(c(-0.1, 1.1)))
})
