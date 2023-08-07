#' Computes autocovariances at lags 0 to p for the ell-th differences of data.
#'
#' @keywords internal
#' @param data       Time series for which we calculate autocovariances.
#' @param ell        Order of differences used.
#' @param p          Maximum lag for the autocovariances
#' @return           Vector of length (p + 1) that consists of empirical
#'                   autocoavariances for the corresponding lag - 1.

emp_acf <- function(data, ell, p) {
  y_diff   <- diff(data, ell)
  len      <- length(y_diff)
  autocovs <- rep(0, p + 1)
  for (k in seq_len(p + 1)) {
    autocovs[k] <- sum(y_diff[k:len] * y_diff[1:(len - k + 1)]) / len
  }
  return(autocovs)
}

#' Computes variance of AR(p) innovation terms eta.
#'
#' @keywords internal
#' @param data      Time series.
#' @param coefs     Estimated coefficients.
#' @param p         AR order of the time series.
#' @return          Variance of the innovation term

variance_eta <- function(data, coefs, p) {
  y_diff <- diff(data)
  len    <- length(y_diff)
  x_diff <- matrix(0, nrow = len - p, ncol = p)
  for (i in seq_len(p)) {
     x_diff[, i] <- y_diff[(p + 1 - i):(len - i)]
  }
  y_diff  <- y_diff[(p + 1):len]
  resid   <- y_diff - as.vector(x_diff %*% coefs)
  var_eta <- mean(resid^2) / 2
  return(var_eta)
}

#' Computes vector of correction terms for second-stage estimator
#' of AR parameters as described in Khismatullina, Vogt (2019).
#'
#' @keywords internal
#' @param coefs     Given (estimated) coefficients of the AR time series.
#' @param var_eta   Variance of the innovation term
#' @param len       Length of the vector of the corrections
#' @return.         Vector of the corrections terms of length len.
corrections <- function(coefs, var_eta, len) {
  p        <- length(coefs)
  c_vec    <- rep(0, len)
  c_vec[1] <- 1
  for (j in 2:len) {
    lags     <- (j - 1):max(j - p, 1)
    c_vec[j] <- sum(coefs[1:length(lags)] * c_vec[lags])
  }
  c_vec <- c_vec * var_eta
  return(c_vec)
}

#' Computes estimator of the AR(p) coefficients by the procedure from
#' Khismatullina and Vogt (2020).
#'
#' @keywords internal
#' @param data      Time series.
#' @param l1,l2     Tuning parameters.
#' @param correct   Vector of the corrections, either zero or calculated by
#'                  the function \code{\link{corrections}}.
#' @param p         AR order of the time series.
#' @return.         Vector of length p  of estimated AR coefficients.
ar_coef <- function(data, l1, l2, correct, p) {
  pos <- 0
  a_mat <- matrix(0, nrow = l2 - l1 + 1, ncol = p)
  for (ell in l1:l2) {
    pos <- pos + 1
    autocovs <- emp_acf(data, ell, p)
    cov_mat  <- matrix(0, ncol = p, nrow = p)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
         cov_mat[i, j] <- autocovs[abs(i - j) + 1]
      }
    }
    if (ell >= p) {
      correct_vec <- correct[(ell - 1 + 1):(ell - p + 1)]
    } else {
      correct_vec <- c(correct[(ell - 1 + 1):1], rep(0, p - ell))
    }
    cov_vec <- autocovs[2:(p + 1)] + correct_vec
    a_mat[pos, ] <- solve(cov_mat) %*% cov_vec
  }
  a_hat <- colMeans(a_mat)
  a_hat <- as.vector(a_hat)
  return(a_hat)
}

#' Computes estimator of the long-run variance of the error terms.
#'
#' @description     A difference based estimator for the coefficients and
#'                  long-run variance in case of a nonparametric regression
#'                  model are AR(p).
#'                  
#'                  Specifically, we assume that we observe \eqn{Y(t)} that satisfy
#'                  the following equation: \deqn{Y(t) = m(t/T) + \epsilon_t.}
#'                  Here, \eqn{m(\cdot)} is an unknown function, and the errors
#'                  \eqn{\epsilon_t} are AR(p) with p known. Specifically, we ler
#'                  \eqn{\{\epsilon_t\}} be a process of the form
#'                  \deqn{\epsilon_t = \sum_{j=1}^p a_j \epsilon_{t-j} + \eta_t,} 
#'                  where \eqn{a_1,a_2,\ldots, a_p} are unknown coefficients and
#'                  \eqn{\eta_t} are i.i.d.\ with \eqn{E[\eta_t] = 0} and
#'                  \eqn{E[\eta_t^2] = \nu^2}.
#'                  
#'                  This function produces an estimator \eqn{\widehat{\sigma}^2}
#'                  of the long-run variance 
#'                  \deqn{\sigma^2 = \sum_{l=-\infty}^{\infty} cov(\epsilon_0,\epsilon_{l})}
#'                  of the error terms, as well as estimators
#'                  \eqn{\widehat{a}_1, \ldots, \widehat{a}_p} of the coefficients
#'                  \eqn{a_1,a_2,\ldots, a_p} and an estimator \eqn{\widehat{\nu}^2} of 
#'                  the innovation variance \eqn{\nu^2}.
#'                  
#'                  The exact estimation procedure as well as description of 
#'                  the tuning parameters needed for this estimation can be found
#'                  in Khismatullina and Vogt (2020).
#'                  
#' @param data      A vector of \eqn{Y(1), Y(2), \ldots, Y(T)}.
#' @param q,r_bar   Tuning parameters.
#' @param p         AR order of the error terms.
#' @return A list with the following elements:
#'     \item{lrv}{Estimator of the long run variance of the error terms
#'     \eqn{\sigma^2}.}
#'     \item{ahat}{Vector of length p of estimated AR coefficients
#'     \eqn{a_1,a_2,\ldots, a_p}.}
#'     \item{vareta}{Estimator of the variance of the innovation term \eqn{\nu^2}.}
#' 
#' @references Khismatullina M., Vogt M. Multiscale inference and long-run
#'             variance estimation in non-parametric regression with
#'             time series errors //Journal of the Royal Statistical Society:
#'             Series B (Statistical Methodology). - 2020.
#'
#' @export
estimate_lrv <- function(data, q, r_bar, p) {
    a_tilde       <- ar_coef(data = data, l1 = q, l2 = q,
                             correct = rep(0, max(p, q) + 1), p = p)
    sig_eta_tilde <- variance_eta(data = data, coefs = a_tilde, p = p)
    correct       <- corrections(coefs = a_tilde, var_eta = sig_eta_tilde,
                                 len = r_bar + 1)
    a_hat         <- ar_coef(data = data, l1 = 1, l2 = r_bar,
                             correct = correct, p = p)
    sig_eta_hat   <- variance_eta(data = data, coefs = a_hat, p = p)
    lrv_hat       <- sig_eta_hat / (1 - sum(a_hat))^2
    return(list(lrv = lrv_hat, ahat = a_hat, vareta = sig_eta_hat))
}
