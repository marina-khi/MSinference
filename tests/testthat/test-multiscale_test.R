test_that("multiscale_test works for single time series", {
  set.seed(101)
  t_len <- 100
  # Simple trend + noise
  trend <- seq(0, 1, length.out = t_len)
  data <- trend + rnorm(t_len, sd = 0.1)
  
  # Use small grid for speed
  u_grid <- seq(0.2, 0.8, by = 0.2)
  h_grid <- c(0.1, 0.2)
  grid <- construct_grid(t_len, u_grid = u_grid, h_grid = h_grid)
  
  result <- multiscale_test(data, sigma = 0.1, n_ts = 1, 
                           grid = grid, sim_runs = 50)
  
  expect_type(result, "list")
  expect_true("statistics" %in% names(result))
  expect_true("quant" %in% names(result))
  expect_true("testing_result" %in% names(result))
  expect_true("test_matrix" %in% names(result))
  
  # Statistics should be finite
  expect_true(is.finite(result$statistics))
  expect_true(is.finite(result$quant))
})

test_that("multiscale_test works for multiple time series", {
  set.seed(202)
  t_len <- 80
  n_ts <- 3
  
  # Generate multiple series
  data <- matrix(rnorm(t_len * n_ts), nrow = t_len, ncol = n_ts)
  
  # Simple grid
  u_grid <- seq(0.25, 0.75, by = 0.25)
  h_grid <- c(0.15, 0.25)
  grid <- construct_grid(t_len, u_grid = u_grid, h_grid = h_grid)
  
  # Compare series 1 vs 2
  ijset <- matrix(c(1, 2), nrow = 1, ncol = 2)
  
  result <- multiscale_test(data, sigma = 1, n_ts = n_ts, 
                           grid = grid, ijset = ijset, 
                           sim_runs = 50, epidem = FALSE)
  
  expect_type(result, "list")
  expect_true("statistics" %in% names(result))
  expect_true("quant" %in% names(result))
  expect_true("stat_pairwise" %in% names(result))
  expect_true(is.finite(result$statistics))
})

test_that("multiscale_test respects alpha parameter", {
  set.seed(303)
  t_len <- 100
  data <- rnorm(t_len)
  
  u_grid <- seq(0.3, 0.7, by = 0.2)
  h_grid <- c(0.15)
  grid <- construct_grid(t_len, u_grid = u_grid, h_grid = h_grid)
  
  result_005 <- multiscale_test(data, sigma = 1, grid = grid, 
                               alpha = 0.05, sim_runs = 100)
  result_001 <- multiscale_test(data, sigma = 1, grid = grid, 
                               alpha = 0.01, sim_runs = 100)
  
  # Quantile at alpha=0.01 should be larger than at alpha=0.05
  expect_true(result_001$quant >= result_005$quant)
})

test_that("multiscale_test validates input dimensions", {
  # Should handle edge cases
  expect_no_error({
    data <- rnorm(50)
    grid <- construct_grid(50, u_grid = c(0.5), h_grid = c(0.15))
    multiscale_test(data, sigma = 1, grid = grid, sim_runs = 10)
  })
})
