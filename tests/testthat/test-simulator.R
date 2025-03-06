test_that("initialize_simulator works correctly", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r, model = "issm")

  expect_s3_class(mod, "tssim.component")
  expect_equal(length(mod$simulated), 50)
  expect_true(all(c("Error") %in% colnames(mod$components)))
})

# Test add_polynomial

test_that("add_polynomial adds correct components", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_polynomial(order = 2, alpha = 0.2, beta = 0.01)

  expect_true(all(c("Level", "Slope") %in% colnames(mod$components)))
  expect_equal(nrow(mod$components), 51)
})

# Test add_seasonal

test_that("add_seasonal adds a seasonal component", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_seasonal(frequency = 12, gamma = 0.5)

  expect_true(any(grepl("Seasonal", colnames(mod$components))))
})

# Test add_arma

test_that("add_arma adds ARMA component", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_arma(order = c(1,1), ar = 0.5, ma = -0.3)

  expect_true("ARMA" %in% colnames(mod$components))
})

# Test add_regressor

test_that("add_regressor adds correct regressor component", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  xreg <- matrix(rnorm(50), ncol = 1)
  mod <- initialize_simulator(r) |> add_regressor(xreg, pars = 0.8)

  expect_true("X" %in% colnames(mod$components))
})

# Test add_transform

test_that("add_transform applies transformation correctly", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_transform(method = "box-cox", lambda = 0.5)

  expect_true(any(mod$table$type == "transform"))
})

# Test add_anomaly

test_that("add_anomaly introduces anomaly correctly", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_anomaly(time = 25, delta = 0.5, ratio = 1.2)
  expect_true(any(grepl("AO|LS|TC", colnames(mod$components))))
})

# Test tsdecompose

test_that("tsdecompose returns correct components", {
  set.seed(123)
  r <- rnorm(50, mean = 0, sd = 5)
  mod <- initialize_simulator(r) |> add_polynomial(order = 2) |> add_seasonal(frequency = 12)
  decomp <- tsdecompose(mod)

  expect_true(ncol(decomp) > 1)
  expect_equal(nrow(decomp), 50)
})

# Test tsensemble

test_that("tsensemble correctly combines multiple simulations", {
  set.seed(123)
  r <- rnorm(100, mean = 0, sd = 5)
  mod1 <- initialize_simulator(r) |> add_polynomial(order = 2)
  mod2 <- initialize_simulator(r) |> add_polynomial(order = 2)
  ens <- mixture_modelspec(mod1, mod2)

  weights <- matrix(c(rep(1, 50), rep(0, 50), rep(0, 50), rep(1, 50)), ncol = 2)
  ens_sim <- tsensemble(ens, weights = weights, difference = FALSE)

  expect_equal(length(ens_sim), 100)
})
