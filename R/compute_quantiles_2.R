#' Computes quantiles of the gaussian multiscale statistics.
#'
#' @description        Quantiles from the gaussian version of the test
#'                     statistics which are used to approximate
#'                     the critical values for the multiscale test.
#' @param t_len        Sample size.
#' @param n_ts         Number of time series analyzed. Default is 1.
#' @param grid         Grid of location-bandwidth points as produced by
#'                     the function \code{\link{construct_grid}} or
#'                     \code{\link{construct_weekly_grid}}, list with
#'                     the elements 'gset', 'bws', 'gtype'. If not provided,
#'                     then the defalt grid is produced and used.
#'                     For the construction of the default grid,
#'                     see \code{\link{construct_grid}}.
#' @param ijset        A matrix of integers. In case of multiple time series,
#'                     we need to know which pairwise comparisons to perform.
#'                     This matrix consists of all pairs of indices \eqn{(i, j)}
#'                     that we want to compare. If not provided, then all
#'                     possible pairwise comparison are performed.
#' @param sigma        Value of \eqn{\sqrt{\sigma^2}}. In case of n_ts = 1,
#'                     \eqn{\sigma^2} denotes the long-run error variance, and
#'                     in case of n_ts > 1, \eqn{\sigma^2} denotes the
#'                     overdispersion parameter.
#'                     If not given, then the default is 1.
#' @param deriv_order  In case of a single time series analysed, this parameter
#'                     denotes the order of the derivative of the trend
#'                     function that is being estimated. Default is 0.
#' @param sim_runs     Number of simulation runs to produce quantiles.
#'                     Default is 1000.
#' @param probs        A numeric vector of probability levels \eqn{(1-\alpha)}
#'                     for which the quantiles are computed.
#'                     Default is \eqn{(0.5, 0.505, 0.51, \ldots, 0.995)}.
#' @param correction   Logical variable, TRUE (by default) if we are using
#'                     \eqn{a_k} and \eqn{b_k}.
#' @param epidem       Logical variable, TRUE if we are using
#'                     dealing with epidemic time trends. Default is FALSE.
#'
#' @return             Matrix with 2 rows where the first row contains
#'                     the vector of probabilities (probs) and the second
#'                     contains corresponding quantiles of the gaussian
#'                     statistics distribution.
#' @export
#'
#' @examples
#' compute_quantiles(100)

compute_quantiles_2 <- function(t_len, n_ts = 1, grid = NULL,
                                ijset = NULL, sigma = 1,
                                deriv_order = 0, sim_runs = 1000,
                                probs = seq(0.5, 0.995, by = 0.005),
                                correction = TRUE, epidem = FALSE) {
  
  numCores  <- round(parallel::detectCores() * .70)
  
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }
  
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }
  
  ijset_cpp               <- as.matrix(ijset)
  ijset_cpp               <- as.vector(ijset_cpp)
  storage.mode(ijset_cpp) <- "integer"
  
  gset                   <- grid$gset
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp)
  storage.mode(gset_cpp) <- "double"
  storage.mode(sigma)    <- "double"
  storage.mode(sigma)    <- "double"

  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
    simulate_gaussian_2(t_len = t_len, n_ts = n_ts, 
                        gset = gset_cpp, ijset = ijset_cpp,
                        sigma = sigma, deriv_order = deriv_order,
                        correction = correction, epidem = epidem)
    # Loop one-by-one using foreach
  } -> phi
  stopCluster(cl)
  
  quant  <- as.vector(quantile(phi, probs = probs))
  quant  <- rbind(probs, quant)
  
  colnames(quant) <- NULL
  rownames(quant) <- NULL
  
  return(list("quant" = quant, "phi" = phi))
}
