#' Multiscale Inference for Nonparametric Time Trend(s)
#'
#'
#' @description This package performs a multiscale analysis of a single nonparametric 
#'   time trends (Khismatullina and Vogt (2020)) or multiple nonparametric 
#'   time trends (Khismatullina and Vogt (2022), Khismatullina and Vogt (2023)).
#'
#'   In case of a single nonparametric regression, the multiscale method to
#'   test qualitative hypotheses about the nonparametric time trend \eqn{m}
#'   in the model \eqn{Y_t = m(t/T) + \epsilon_t} with time series errors
#'   \eqn{\epsilon_t} is provided. The method was first proposed in
#'   Khismatullina and Vogt (2020). It allows to test for shape properties
#'   (areas of monotonic decrease or increase) of the trend \eqn{m}.
#'
#'   This method require an estimator of the long-run error variance
#'   \eqn{\sigma^2 = \sum_{l=-\infty}^{\infty} Cov(\epsilon_0, \epsilon_l)}.
#'   Hence, the package also provides the difference-based
#'   estimator for the case that the errors belong to the class of
#'   \eqn{AR(\infty)} processes. The estimator was also proposed in
#'   Khismatullina and Vogt (2020).
#'
#'   In case of multiple nonparametric regressions, we provide
#'   the multiscale method to test qualitative hypotheses about
#'   the nonparametric time trends in the context of epidemic modelling.
#'   Specifically, we assume that the we observe a sample of the count data
#'   \eqn{\{\mathcal{X}_i = \{ X_{it}: 1 \le 1 \le T \}\}}, where \eqn{X_{it}}
#'   are quasi-Poisson distributed with time-varying intensity parameter
#'   \eqn{\lambda_i(t/T)}. The multiscale method allows to test whether
#'   intensity parameters are different or not, and if they are, it detects
#'   with a pre-specified significance level the regions where these differences
#'   most probably occur. The method was introduced in
#'   Khismatullina and Vogt (2023) and can be used for comparing the rates of
#'   infection of COVID-19 across countries.
#'
#' @name MSinference-package
#' @aliases MSinference
#' @references 
#' \insertRef{KhismatullinaVogt2020}{MSinference}
#'
#' \insertRef{KhismatullinaVogt2022}{MSinference}            
#'             
#' \insertRef{KhismatullinaVogt2023}{MSinference}            
NULL
