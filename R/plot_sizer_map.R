#' Plots SiZer map from the test results of the multiscale testing procedure.
#'
#' @param u_grid       Vector of location points in the unit interval \eqn{[0,1]}.
#' @param h_grid       Vector of bandwidths from \eqn{(0,0.5)}.
#' @param test_results Matrix of test results created by
#'                     \code{\link{multiscale_test}}.
#' @param plot_title   Title of the plot. Default is NA and no title is written.
#' @param greyscale    Whether SiZer map is plotted in grey scale.
#'                     Default is FALSE.
#' @param ...          Any further options to be passed to the image function.
#'
#' @export

plot_sizer_map <- function(u_grid, h_grid, test_results, plot_title = NA,
                           greyscale = FALSE, ...) {

  if (greyscale) {
    col_vec <- c("#F7F7F7", "#969696", "#525252", "#636363")
  } else {
    col_vec <- c("red", "purple", "blue", "gray")
  }
  temp    <- sort(unique(as.vector(test_results))) + 2
  temp    <- seq(min(temp), max(temp), by = 1)
  col_vec <- col_vec[temp]

  image(x = u_grid, y = log(h_grid, 10), z = t(test_results), col = col_vec,
        xlab = "", ylab = expression(log[10](h)), main = plot_title,
        xaxt = "n", mgp = c(1, 0.5, 0), ...)
}
