test_that("emp_acf computes autocovariances correctly", {
  # Simple test with known series
  set.seed(42)
  data <- rnorm(100)
  ell <- 1
  p <- 5
  
  result <- MSinference:::emp_acf(data, ell, p)
  
  expect_length(result, p + 1)
  expect_type(result, "double")
  
  # Autocovariance at lag 0 should be non-negative
  expect_true(result[1] >= 0)
})

test_that("variance_eta computes positive variance", {
  set.seed(123)
  # AR(1) process: y_t = 0.5 * y_{t-1} + e_t
  n <- 200
  data <- arima.sim(list(ar = 0.5), n = n)
  p <- 1
  coefs <- c(0.5)
  
  var_result <- MSinference:::variance_eta(data, coefs, p)
  
  expect_type(var_result, "double")
  expect_length(var_result, 1)
  expect_true(var_result > 0)
})

test_that("corrections returns vector of correct length", {
  coefs <- c(0.5, 0.3)
  var_eta <- 1.0
  len <- 10
  
  c_vec <- MSinference:::corrections(coefs, var_eta, len)
  
  expect_length(c_vec, len)
  expect_type(c_vec, "double")
  expect_equal(c_vec[1], var_eta)  # First element should be var_eta
})

test_that("ar_coef estimates coefficients for AR process", {
  set.seed(456)
  # Generate AR(2) data
  n <- 300
  data <- arima.sim(list(ar = c(0.6, -0.3)), n = n)
  
  p <- 2
  l1 <- 1
  l2 <- 3
  correct <- rep(0, p)
  
  coefs <- MSinference:::ar_coef(data, l1, l2, correct, p)
  
  expect_length(coefs, p)
  expect_type(coefs, "double")
  # Coefficients should be finite
  expect_true(all(is.finite(coefs)))
})

test_that("estimate_lrv returns positive values", {
  set.seed(789)
  data <- cumsum(rnorm(200))  # Random walk with drift
  
  result <- estimate_lrv(data)
  
  expect_type(result, "list")
  expect_true("sigma" %in% names(result))
  expect_true(result$sigma > 0)
  expect_true(is.finite(result$sigma))
})
