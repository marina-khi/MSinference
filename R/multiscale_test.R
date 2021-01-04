#' Carries out the multiscale test given that the values the estimatates of
#' long-run variance have already been computed.
#'
#' @param n_ts          Number of time series analysed. Default is 1.
#' @param data          Vector (in case of n_ts = 1) or matrix (in case of
#'                      n_ts > 1) that contains (a number of) time series
#'                      that needs to be analyzed. In the latter case,
#'                      each column of the matrix must contain one time series.
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
#' @param alpha        Significance level. Default is \eqn{0.05}.
#' @param sim_runs     Number of simulation runs to produce quantiles.
#'                     Default is 1000.
#' @param deriv_order  In case of a single time series, this denotes the order of
#'                     the derivative of the trend that we estimate.
#'                     Default is 0.
#' @param correction   Logical variable, TRUE (by default) is we are using
#'                     \eqn{a_k} and \eqn{b_k}.
#' @export
#'
#' @return In case of n_ts = 1, the function returns a list
#' with the following elements:
#'    \item{quant}{Quantile that was used for testing calculated from
#'                the gaussian distribution.}
#'    \item{statistics}{Value of the multiscale statistics.}
#'    \item{test_matrix}{Matrix of the test results for the multiscale test
#'                       defined in Khismatullina and Vogt (2019).
#'                       The matrix is coded as follows:
#'                       \itemize{
#'                       \item test_matrix[i,j] = -1: test rejects the null for the
#'                                j-th location \eqn{u} and the i-th bandwidth \eqn{h} and
#'                                indicates a decrease in the trend;
#'                       \item test_matrix[i,j] = 0: test does not reject the null
#'                                for the j-th location \eqn{u} and the i-th
#'                                bandwidth \eqn{h};
#'                       \item test_matrix[i,j] = 1:  test rejects the null for the
#'                                j-th location \eqn{u} and the i-th bandwidth \eqn{h} and
#'                                indicates an increase in the trend;
#'                       \item test_matrix[i,j] = 2: no test is carried out at j-th
#'                                location \eqn{u} and i-th bandwidth \eqn{h} (because
#'                                the point \eqn{(u, h)} is excluded from the grid
#'                                as specified by the 'deletions' option
#'                                in the function \code{\link{construct_grid}})}.
#'                       }
#'    \item{gset_with_vals}{A matrix that contains the values of the normalised 
#'                         kernel averages and test results for each pair
#'                         of location-bandwidth
#'                         with the corresponding location and bandwidth.}
#' In case of n_ts > 1, the function returns a list
#' with the following elements:
#'    \item{quant}{Quantile that was used for testing calculated from
#'                the gaussian distribution.stat}{Value of the multiscale statistics.}
#'    \item{statistics}{Value of the multiscale statistics.}
#'    \item{stat_pairwise}{Matrix of the values of the pairwise statistics.}
#'    \item{ijset}{The matrix that  consists of all pairs of indices
#'                 \eqn{(i, j)} that we compared. The order of these
#'                 pairs corresponds to the order in the list
#'                 gset_with_vals.}
#'    \item{gset_with_vals}{A list of matrices, each matrix corresponding to a 
#'                         specific pairwise comparison. The order of the list 
#'                         is determined by ijset. Each matrix contains
#'                         the values of the normalisedkernel averages
#'                         for each pair of location-bandwidth
#'                         with the corresponding location and bandwidth.}
multiscale_test <- function(data, sigma, n_ts = 1, grid = NULL,
                            ijset = NULL, alpha = 0.05, sim_runs = 1000,
                            deriv_order = 0, correction = TRUE) {

  if (n_ts == 1) {
    t_len <- length(data)
  } else {
    t_len <- nrow(data)
  }

  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }

  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }


  #These full grids we need for plotting the SiZer maps
  gset_full    <- grid$gset_full
  u_grid_full  <- unique(gset_full[, 1])
  pos_full     <- grid$pos_full

  test_res     <- rep(2, length(pos_full))

  # Select (1-alpha) quantile of the multiscale statistic under the null
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sigma = sigma,
                                 sim_runs = sim_runs,
                                 deriv_order = deriv_order,
                                 correction = correction)

  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])

  if (sum(probs == (1 - alpha)) == 0)
    pos <- which.min(abs(probs - (1 - alpha)))
  if (sum(probs == (1 - alpha)) != 0)
    pos <- which.max(probs == (1 - alpha))

  quant <- quant[pos]

  psi   <- compute_statistics(data = data, sigma = sigma, n_ts = n_ts,
                              grid = grid, deriv_order = deriv_order)
  stat  <- psi$stat

  if (n_ts == 1) {
    gset_with_vals             <- psi$gset_with_vals
    test_results               <- (gset_with_vals$vals_cor > quant) * sign(gset_with_vals$vals)
    test_res[!is.na(pos_full)] <- test_results
    test                       <- matrix(test_res, ncol = length(u_grid_full),
                                         byrow = TRUE)

    if (stat > quant) {
      cat("For the given time series we reject H_0 with probability",
          alpha, ".\n Psihat_statistic = ", stat,
          ".\n Gaussian quantile value = ", quant, "\n")
    } else {
      cat("For the given time series we fail to reject H_0 with probability",
          alpha, ".\n Psihat_statistic = ", stat,
          ".\n Gaussian quantile value = ", quant, "\n")
    }

    gset_with_vals$test <- test_results

    return(list(quant = quant, stat = psi$stat, test_matrix = test,
                gset_with_vals = gset_with_vals))

  } else {
    if (stat > quant) {
      cat("We reject H_0 with probability", alpha, ". Psihat_statistic = ",
          stat, ".\n Number of pairwise rejections = ",
          sum(psi$stat_pairwise > quant, na.rm = TRUE), "out of ", nrow(ijset),
          ".\n Gaussian quantile value = ", quant, "\n")
    } else {
      cat("We fail to reject H_0 with probability", alpha,
          ".\n Psihat_statistic = ", stat,
          ".\n Gaussian quantile value = ", quant, "\n")
    }

    gset_with_values <- psi$gset_with_values

    for (i in seq_len(nrow(ijset))) {
      if (correction) {
        test_results               <- gset_with_values[[i]]$vals_cor > quant
      } else {
        test_results               <- gset_with_values[[i]]$vals > quant
      }
      gset_with_values[[i]]$test <- test_results
    }

    return(list(quant = quant, stat = stat, stat_pairwise = psi$stat_pairwise,
                ijset = ijset, gset_with_values = gset_with_values))
  }
}
