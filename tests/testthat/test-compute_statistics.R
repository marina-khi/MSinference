test_that("compute_statistics works for single time series", {
  set.seed(111)
  t_len <- 100
  data <- cumsum(rnorm(t_len))
  
  u_grid <- seq(0.2, 0.8, by = 0.2)
  h_grid <- c(0.1, 0.2)
  grid <- construct_grid(t_len, u_grid = u_grid, h_grid = h_grid)
  gset <- cbind(grid$gset[, 1], grid$gset[, 2])
  gset_vec <- c(gset[, 1], gset[, 2])
  
  sigma <- 1.0
  deriv_order <- 0
  
  result <- compute_statistics(t_len, data, gset_vec, sigma, deriv_order)
  
  expect_type(result, "list")
  expect_true("stat" %in% names(result))
  expect_true("vals" %in% names(result))
  expect_true("vals_cor" %in% names(result))
  
  expect_true(is.finite(result$stat))
  expect_length(result$vals, nrow(gset))
  expect_length(result$vals_cor, nrow(gset))
})

test_that("compute_statistics handles different derivative orders", {
  set.seed(222)
  t_len <- 120
  data <- seq(0, 1, length.out = t_len) + rnorm(t_len, sd = 0.1)
  
  grid <- construct_grid(t_len, u_grid = c(0.5), h_grid = c(0.2))
  gset_vec <- c(grid$gset[, 1], grid$gset[, 2])
  
  result0 <- compute_statistics(t_len, data, gset_vec, 
                                sigma = 1, deriv_order = 0)
  result1 <- compute_statistics(t_len, data, gset_vec, 
                                sigma = 1, deriv_order = 1)
  
  expect_true(is.finite(result0$stat))
  expect_true(is.finite(result1$stat))
  # Different derivative orders should give different results
  expect_false(isTRUE(all.equal(result0$vals, result1$vals)))
})

test_that("compute_quantiles produces reasonable quantiles", {
  set.seed(333)
  t_len <- 100
  sim_runs <- 100
  
  grid <- construct_grid(t_len, u_grid = c(0.5), h_grid = c(0.15))
  gset_vec <- c(grid$gset[, 1], grid$gset[, 2])
  
  ijset <- c(1, 1)
  alpha <- 0.05
  
  result <- compute_quantiles(t_len, n_ts = 1, sim_runs = sim_runs,
                             gset = gset_vec, ijset = ijset,
                             sigma = 1, deriv_order = 0,
                             alpha = alpha)
  
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result > 0)
  expect_true(is.finite(result))
})

test_that("compute_minimal_intervals identifies intervals correctly", {
  # Create test matrix with known structure
  test_matrix <- matrix(0, nrow = 3, ncol = 5)
  test_matrix[2, 2:4] <- 1  # Rejection region in middle
  
  gset <- expand.grid(u = seq(0.2, 1, by = 0.2), h = c(0.1, 0.15, 0.2))
  
  result <- compute_minimal_intervals(test_matrix, gset$u, gset$h)
  
  expect_type(result, "list")
  # Result should contain interval information
  expect_true(length(result) > 0)
})
