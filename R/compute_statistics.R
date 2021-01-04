#' Calculates the value of the test statistics both for single time series
#' analysis and multiple time series analysis.
#'
#' @param n_ts         Number of time series analysed. Default is 1.
#' @param data         Vector (in case of n_ts = 1) or matrix (in case of
#'                     n_ts > 1) that contains (a number of) time series
#'                     that needs to be analyzed. In the latter case,
#'                     each column of the matrix must contain one time series.
#' @param sigma        The estimator of the square root of the long-run
#'                     variance \eqn{\sigma} in case of n_ts = 1,
#'                     or the estimator of the overdispersion parameter
#'                     \eqn{\sigma} in case of n_ts > 1.
#' @param grid         Grid of location-bandwidth points as produced by
#'                     the functions \code{\link{construct_grid}} or
#'                     \code{\link{construct_weekly_grid}}, it is a list with
#'                     the elements 'gset', 'bws', 'gtype'. If not provided,
#'                     then the defalt grid is used.
#'                     For the construction of the default grid,
#'                     see \code{\link{construct_grid}}.
#' @param ijset        In case of multiple time series (n_ts > 1),
#'                     we need to know which pairs of time series to compare.
#'                     This matrix consists of all pairs of indices \eqn{(i, j)}
#'                     that we want to compare. If not provided, then all
#'                     possible pairwise comparison are performed.
#' @param deriv_order  In case of a single time series, this denotes the order of
#'                     the derivative of the trend that we estimate.
#'                     Default is 0.
#' @export
#'
#' @return In case of n_ts = 1, the function returns a list
#' with the following elements:
#' \item{stat}{Value of the multiscale statistics.}
#' \item{gset_with_vals}{A matrix that contains the values of the normalised 
#'                         kernel averages for each pair of location-bandwidth
#'                         with the corresponding location and bandwidth.}
#'                         
#' In case of n_ts > 1, the function returns a list
#' with the following elements:
#' \item{stat}{Value of the multiscale statistics.}
#' \item{stat_pairwise}{Matrix of the values of the pairwise statistics.}
#' \item{ijset}{The matrix that  consists of all pairs of indices
#'                         \eqn{(i, j)} that we compared. The order of these
#'                         pairs corresponds to the order in the list
#'                         gset_with_vals.}
#' \item{gset_with_vals}{A list of matrices, each matrix corresponding to a 
#'                         specific pairwise comparison. The order of the list 
#'                         is determined by ijset. Each matrix contains
#'                         the values of the normalisedkernel averages
#'                         for each pair of location-bandwidth
#'                         with the corresponding location and bandwidth.}
compute_statistics <- function(data, sigma, n_ts = 1, grid = NULL,
                            ijset = NULL, deriv_order = 0) {

  if (n_ts == 1) {
    t_len <- length(data)
  } else {
    t_len <- nrow(data)
  }

  t_len <- as.integer(t_len)

  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }

  # Compute test results
  gset                   <- grid$gset
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp)
  storage.mode(gset_cpp) <- "double"

  if (n_ts == 1) {
    psi   <- compute_single_statistics(t_len = t_len, data = data,
                                       gset = gset_cpp, sigma = sigma,
                                       deriv_order = deriv_order)

    gset$vals     <- as.vector(psi$vals)
    gset$vals_cor <- as.vector(psi$vals_cor)

    return(list(stat = psi$stat, gset_with_vals = gset))

  } else {
    if (is.null(ijset)) {
      ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
      ijset <- ijset[ijset$i < ijset$j, ]
    }

    ijset_cpp               <- as.matrix(ijset)
    ijset_cpp               <- as.vector(ijset_cpp)
    storage.mode(ijset_cpp) <- "integer"

    psi_ij <- compute_multiple_statistics(t_len = t_len, n_ts = n_ts,
                                          data = data, gset = gset_cpp,
                                          ijset = ijset_cpp,
                                          sigma = sigma)
    stat_max <- max(psi_ij$stat, na.rm = TRUE)
    gset_with_values <- list()

    for (i in seq_len(nrow(ijset))) {
      gset$vals              <- psi_ij$vals_matrix[, i]
      gset$vals_cor          <- psi_ij$vals_cor_matrix[, i]
      gset_with_values[[i]]  <- gset
    }

    return(list(stat = stat_max, stat_pairwise = psi_ij$stat,
                ijset = ijset, gset_with_values = gset_with_values))
  }
}
