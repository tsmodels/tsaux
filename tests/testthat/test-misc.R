test_that("seasonality test (weeks)", {
  set.seed(123)
  r <- rnorm(52*4, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_seasonal(frequency = 52, gamma = 0.5, init_scale = 5)
  d <- future_dates(as.Date("2000-01-01"), n = 52*4, frequency = "weeks")
  x <- xts(mod$simulated, d)
  expect_true(seasonality_test(x))
})


test_that("seasonality test (months)/1", {
  set.seed(123)
  r <- rnorm(12*4, mean = 0, sd = 1)
  mod <- initialize_simulator(r) |> add_polynomial(alpha = 0.9) |> add_seasonal(frequency = 12, gamma = 0.5, init_scale = 3)
  d <- future_dates(as.Date("2000-01-01"), n = 12*4, frequency = "months")
  x <- xts(mod$simulated, d)
  expect_false(seasonality_test(x))
})


test_that("seasonality test (months)/2", {
  set.seed(123)
  r <- rnorm(12*5, mean = 0, sd = 1)
  mod <- initialize_simulator(r) |> add_polynomial(alpha = 0.9) |> add_seasonal(frequency = 12, gamma = 0.5, init_scale = 3)
  d <- future_dates(as.Date("2000-01-01"), n = 12*5, frequency = "months")
  x <- xts(mod$simulated, d)
  expect_true(seasonality_test(x))
})


test_that("seasonality test (hours)", {
  set.seed(123)
  r <- rnorm(24*7*5, mean = 0, sd = 1)
  mod <- initialize_simulator(r) |> add_polynomial(alpha = 0.9) |> add_seasonal(frequency = 24, gamma = 0.5, init_scale = 3)
  d <- future_dates(as.Date("2000-01-01"), n = 24*7*5, frequency = "hours")
  x <- xts(mod$simulated, d)
  expect_true(seasonality_test(x))
})
