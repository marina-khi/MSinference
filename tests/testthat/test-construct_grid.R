test_that("construct_grid creates default grid correctly", {
  t_len <- 100
  grid <- construct_grid(t_len)
  
  # Check return structure
  expect_type(grid, "list")
  expect_named(grid, c("gset", "bws", "lens", "gtype", "gset_full", "pos_full"))
  
  # Check grid type
  expect_equal(grid$gtype, "default")
  
  # Check dimensions
  expect_true(nrow(grid$gset) > 0)
  expect_equal(ncol(grid$gset), 2)
  expect_equal(ncol(grid$gset_full), 2)
  
  # Check bandwidth constraints (0 < h < 0.5)
  expect_true(all(grid$gset[, 2] > 0))
  expect_true(all(grid$gset[, 2] < 0.5))
  
  # Check location constraints (0 < u <= 1)
  expect_true(all(grid$gset[, 1] > 0))
  expect_true(all(grid$gset[, 1] <= 1))
})

test_that("construct_grid accepts custom u_grid and h_grid", {
  t_len <- 100
  u_custom <- seq(0.1, 0.9, by = 0.1)
  h_custom <- c(0.1, 0.2, 0.3)
  
  grid <- construct_grid(t_len, u_grid = u_custom, h_grid = h_custom)
  
  expect_equal(grid$gtype, "non-default")
  expect_equal(nrow(grid$gset), length(u_custom) * length(h_custom))
  expect_true(all(grid$bws %in% h_custom))
})

test_that("construct_grid handles deletions correctly", {
  t_len <- 100
  u_custom <- c(0.25, 0.5, 0.75)
  h_custom <- c(0.1, 0.2)
  n_points <- length(u_custom) * length(h_custom)
  
  # Delete first and last points
  deletions <- rep(TRUE, n_points)
  deletions[c(1, n_points)] <- FALSE
  
  grid <- construct_grid(t_len, u_grid = u_custom, h_grid = h_custom, 
                        deletions = deletions)
  
  expect_equal(nrow(grid$gset), n_points - 2)
  expect_equal(nrow(grid$gset_full), n_points)
  expect_equal(sum(grid$pos_full), n_points - 2)
})

test_that("construct_grid validates input dimensions", {
  # Should work with reasonable inputs
  expect_no_error(construct_grid(50))
  expect_no_error(construct_grid(1000))
})
