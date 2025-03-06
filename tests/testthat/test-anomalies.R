test_that("auto_regressors returns a list with correct elements", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  f <- auto_regressors(y, frequency = 1, sampling = "days", lambda = NULL, h = 20, method = "full")

  expect_type(f, "list")
  expect_named(f, c("xreg", "init"))
  expect_true(is.null(f$xreg) || is.xts(f$xreg))  # xreg can be NULL or xts
  expect_true(is.null(f$init) || is.numeric(f$init))  # init can be NULL or numeric vector
})

test_that("auto_regressors detects anomalies", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  # Add an artificial anomaly
  y[20] <- y[20] * 1.5

  f <- auto_regressors(y, frequency = 1, sampling = "days", lambda = NULL, h = 0, method = "full", maxit.oloop = 10, maxit.iloop = 40)

  expect_true(!is.null(f$xreg))  # Should detect anomalies and not return NULL
  expect_true(ncol(f$xreg) > 0)  # At least one detected anomaly
})

test_that("auto_clean removes anomalies correctly", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  # Introduce a strong additive outlier
  y[30] <- y[30] * 2

  # Clean anomalies
  y_cleaned <- auto_clean(y, frequency = 1, lambda = NULL, method = "full", maxit.oloop = 10, maxit.iloop = 40)

  expect_true(is.xts(y_cleaned))
  expect_equal(length(y_cleaned), length(y))
  expect_true(abs(as.numeric(y_cleaned[30]) - as.numeric(y[29])) < abs(as.numeric(y[30]) - as.numeric(y[29]))) # Check if anomaly is reduced
})


test_that("additive_outlier correctly adds anomalies", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  y_outlier <- additive_outlier(y, time = 20, parameter = 0.2)

  expect_true(as.numeric(y_outlier[20]) > as.numeric(y[20]))
  expect_equal(length(y_outlier), length(y))
})

test_that("temporary_change applies decay correctly", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  y_tc <- temporary_change(y, time = 30, parameter = 0.2, alpha = 0.6)

  expect_true(as.numeric(y_tc[30]) > as.numeric(y[30]))  # Change should be present
  expect_true(as.numeric(y_tc[35]) > as.numeric(y[35]))  # Effect should decay over time
  expect_equal(length(y_tc), length(y))
})

test_that("level_shift applies shifts correctly", {
  library(xts)
  set.seed(42)
  y <- xts(cumprod(c(100, (1 + rnorm(100, 0.01, 0.005)))), as.Date(1:101, origin = "2000-01-01"))

  y_ls <- level_shift(y, time = 50, parameter = 0.3)

  expect_true(as.numeric(y_ls[50]) > as.numeric(y[50]))  # Shift should be present
  expect_true(as.numeric(y_ls[60]) > as.numeric(y[60]))  # Shift should persist
  expect_equal(length(y_ls), length(y))
})
