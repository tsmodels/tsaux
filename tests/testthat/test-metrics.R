actual <- c(100, 200, 300, 400, 500)
predicted <- c(110, 190, 290, 410, 490)
weights <- c(0.1, 0.2, 0.2, 0.3, 0.2)

set.seed(42)
distribution <- matrix(rnorm(500, mean = rep(actual, each = 100), sd = 15), nrow = 100, ncol = length(actual))

test_that("MAPE is computed correctly", {
  expected_mape <- mean(abs((actual - predicted) / actual))
  expect_equal(mape(actual, predicted), expected_mape, tolerance = 1e-6)
})

test_that("RMAPE does not return errors", {
  expect_silent(rmape(actual, predicted))
})

test_that("SMAPE is computed correctly", {
  expected_smape <- mean(2 * abs(actual - predicted) / (abs(actual) + abs(predicted)))
  expect_equal(smape(actual, predicted), expected_smape, tolerance = 1e-6)
})

test_that("MASE handles basic case correctly", {
  original_series <- c(90, 180, 270, 360, 450, 540)
  frequency <- 1
  expect_silent(mase(actual, predicted, original_series, frequency))
})

test_that("MSLRE is computed correctly", {
  valid <- actual > 0 & predicted > 0
  expected_mslre <- mean((log(1 + actual[valid]) - log(1 + predicted[valid]))^2)
  expect_equal(mslre(actual, predicted), expected_mslre, tolerance = 1e-6)
})

test_that("Bias is computed correctly", {
  expected_bias <- mean(predicted - actual)
  expect_equal(bias(actual, predicted), expected_bias, tolerance = 1e-6)
})

test_that("WAPE handles weights correctly", {
  expected_wape <- sum(weights * abs(predicted - actual) / actual)
  expect_equal(wape(actual, predicted, weights), expected_wape, tolerance = 1e-6)
})

test_that("WSLRE handles weights correctly", {
  valid <- actual > 0 & predicted > 0
  expected_wslre <- sum(weights * (log(predicted[valid] / actual[valid]))^2)
  expect_equal(wslre(actual, predicted, weights), expected_wslre, tolerance = 1e-6)
})

test_that("WSE handles weights correctly", {
  expected_wse <- sum(weights * (predicted / actual)^2)
  expect_equal(wse(actual, predicted, weights), expected_wse, tolerance = 1e-6)
})

test_that("Pinball Loss handles sample x horizon distribution correctly", {
  tau <- 0.1
  expected_loss <- mean(ifelse(actual >= apply(distribution, 2, quantile, probs = tau),
                               tau * (actual - apply(distribution, 2, quantile, probs = tau)),
                               (1 - tau) * (apply(distribution, 2, quantile, probs = tau) - actual)))
  expect_equal(pinball(actual, distribution, tau), expected_loss, tolerance = 1e-6)
})

test_that("CRPS handles sample x horizon distribution correctly", {
  # Check CRPS runs correctly and does not produce errors
  expect_silent(crps(actual, distribution))
})

test_that("MIS does not return errors", {
  lower <- actual - 20
  upper <- actual + 20
  alpha <- 0.05
  expect_silent(mis(actual, lower, upper, alpha))
})

test_that("MSIS handles sample x horizon distribution correctly", {
  lower <- actual - 20
  upper <- actual + 20
  original_series <- c(90, 180, 270, 360, 450, 540)
  frequency <- 1
  alpha <- 0.05
  expect_silent(msis(actual, lower, upper, original_series, frequency, alpha))
})
