#' Calculates different information criterions for a single time series
#' or multiple time series with AR(\eqn{p}) errors
#' based on the long-run variance estimator(s) for a range of tuning
#' parameters and different orders \eqn{p}.
#'
#' @export
#'
#' @description This function fits AR(1), ... AR(9) models for all
#'              given time series and calculates different information
#'              criterions (FPE, AIC, AICC, SIC, HQ) for each of these fits.
#'              The result is the best fit in terms of minimizing
#'              the infromation criteria.
#' @param data  One or a number of time series in a matrix. Column names
#'              of the matrix should be reasonable
#' @param q     A vector of integers that consisits of different tuning
#'              parameters to analyse. If not supplied, q is taken to be
#'              \eqn{[2\log{T}]:([2\sqrt{T}] + 1)}.
#' @param r     A vector of integers that consisits of different tuning
#'              parameters r_bar for \code{\link{estimate_lrv}}.
#'              If not supplied, \eqn{r = 5, \ldots, 15}.
#' @return      A list with a number of elements:
#' \item{orders}{A vector of chosen orders of length equal to the number
#'          of time series.
#'          For each time series the order is calculated as
#'          \eqn{\max(which.min(FPE), ... which.min(HQ))}}
#' \item{...}{Matrices with the orders that were selected
#'            (among \eqn{1, \ldots, 9}) for each information criterion.
#'            One matrix for each time series.}
select_order <- function(data, q = NULL, r = 5:15) {
  t_len <- nrow(data)
  n_ts <- ncol(data)

  if (is.null(q)) {
    q <- seq(floor(2 * log(t_len)), ceiling(2 * sqrt(t_len)), by = 1)
  }

  result_list   <- list()
  order_results <- c()
  if (n_ts == 1) {
    list_names <- c("the only time series")
  } else {
    list_names <- colnames(data)
  }

  for (j in 1:n_ts) {
    criterion_matrix <- expand.grid(q = q, r = r)

    criterion_matrix$fpe  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$aic  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$aicc <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$sic  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$hq   <- numeric(length = nrow(criterion_matrix))

    for (i in seq_len(nrow(criterion_matrix))) {
      fpe <- c()
      aic <- c()
      aicc <- c()
      sic <- c()
      hq <- c()

      different_orders <- (1:9)

      for (order in different_orders) {
        ar_struc <- estimate_lrv(data = data[, j], q = criterion_matrix$q[[i]],
                                 r_bar = criterion_matrix$r[[i]], p = order)
        sigma_eta_hat <- sqrt(ar_struc$vareta)
        fpe <- c(fpe, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
        aic <- c(aic, t_len * log(sigma_eta_hat^2) + 2 * order)
        aicc <- c(aicc, t_len * log(sigma_eta_hat^2) +
                      t_len * (1 + order / t_len) / (1 - (order + 2) / t_len))
        sic <- c(sic, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
        hq <- c(hq, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
      }
      criterion_matrix$fpe[[i]]  <- which.min(fpe)
      criterion_matrix$aic[[i]]  <- which.min(aic)
      criterion_matrix$aicc[[i]] <- which.min(aicc)
      criterion_matrix$sic[[i]]  <- which.min(sic)
      criterion_matrix$hq[[i]]   <- which.min(hq)
    }
    maxim <- max(criterion_matrix[, 3:7])
    order_results <- c(order_results, maxim)
    cat("For ", list_names[j], " the results are as follows: ",
        max(criterion_matrix$fpe), " ", max(criterion_matrix$aic), " ",
        max(criterion_matrix$aicc), " ", max(criterion_matrix$sic), " ",
        max(criterion_matrix$hq), " \n")
    result_list[[list_names[j]]] <- criterion_matrix
  }
  result_list[["orders"]] <- order_results
  return(result_list)
}
