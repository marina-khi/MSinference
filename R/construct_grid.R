#' Computes the location-bandwidth grid for the multiscale test.
#'
#' @param t          Sample size.
#' @param u_grid     Vector of location points in the unit interval
#'                   \eqn{[0,1]}. If NULL, a default grid is used.
#' @param h_grid     Vector of bandwidths, each bandwidth is supposed to lie
#'                   in \eqn{(0, 0.5)}. If NULL, a default grid is used.
#' @param deletions  Logical vector of the length len(u.grid) * len(h.grid).
#'                   Each element is either TRUE, which means that
#'                   the corresponding location-bandwidth point \eqn{(u, h)}
#'                   is NOT deleted from the grid, or FALSE, which means that
#'                   the corresponding location-bandwidth point \eqn{(u, h)}
#'                   IS deleted from the grid. Default is NULL
#'                   in which case nothing is deleted.
#'                   See vignette for the use.
#' @return A list with the following elements:
#'    \item{gset}{Matrix of location-bandwidth points \eqn{(u, h)}
#'                   that remains after deletions, the i-th row gset[i,]
#'                   corresponds to the i-th point \eqn{(u,h)}.}
#'     \item{bws}{Vector of bandwidths (after deletions).}
#'     \item{lens}{Vector of length = length(bws), lens[i] gives
#'                   the number of locations in the grid for
#'                   the i-th bandwidth level.}
#'     \item{gtype}{Type of grid that is used, either 'default' or
#'                   'non-default'.}
#'     \item{gset_full}{Matrix of all location-bandwidth pairs \eqn{(u, h)}
#'                   including deleted ones.}
#'     \item{pos_full}{Logical vector indicating which points \eqn{(u, h)}
#'                   have been deleted.}
#'
#' @export
#'
#' @examples
#' construct_grid(100)
#' construct_grid(100, u_grid = seq(from = 0.05, to = 1, by = 0.05),
#'                h_grid = c(0.1, 0.2, 0.3, 0.4))


construct_grid <- function(t, u_grid = NULL, h_grid = NULL, deletions = NULL) {

  grid_type <- "non-default"
  if (is.null(u_grid) & is.null(h_grid)) {
    grid_type <- "default"
    t_len  <- t
    u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
    h_grid <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
    h_grid <- h_grid[h_grid > log(t_len) / t_len]
  }

  gset      <- expand.grid(u = u_grid, h = h_grid)
  gset_full <- gset
  n         <- dim(gset)[1]
  pos_full  <- rep(TRUE, n)

  if (!is.null(deletions))
    pos_full <- deletions

  gset <- gset[pos_full, ]
  bws  <- unique(gset[, 2])
  lengths_u <- rep(NA, length(bws))
  for (i in seq_len(length(bws)))
     lengths_u[i] <- sum(gset[, 2] == bws[i])

  return(list(gset = gset, bws = bws, lens = lengths_u, gtype = grid_type,
              gset_full = gset_full, pos_full = pos_full))
}


#' Computes the location-bandwidth weekly grid for the multiscale test.
#'
#' @param t            Sample size.
#' @param min_len      Minimal length of the interval considered. The grid then
#'                     consists of intervals with lengths min_len,
#'                     2 * min_len, 3 * min_len, ... Default is 7, i.e. a week.
#' @param nmbr_of_wks  Number that determines the longest intervals in the grid:
#'                     the length of this interval is calculated then as 
#'                     min_len * nmbr_of_wks. Default is 4.
#' @return A list with the following elements:
#'    \item{gset}{Matrix of location-bandwidth points \eqn{(u, h)}
#'                the i-th row gset[i,] corresponds to the i-th point
#'                \eqn{(u,h)}.}
#'     \item{bws}{Vector of bandwidths.}
#'     \item{lens}{Vector of length = length(bws), lens[i] gives
#'                   the number of locations in the grid for
#'                   the i-th bandwidth level.}
#'     \item{gtype}{Type of grid that is used, always 'default'.}
#'     \item{gset_full}{Matrix of all location-bandwidth pairs \eqn{(u, h)}.}
#'
#' @export
#'
#' @examples
#' construct_weekly_grid(100)
#' construct_weekly_grid(100, min_len = 7, nmbr_of_wks = 2)

construct_weekly_grid <- function(t, min_len = 7, nmbr_of_wks = 4) {
  mat_base <- numeric(0)
  step     <- floor(min_len / 2)
  nbr      <- floor(t / min_len)
  for (k in seq_len(nmbr_of_wks)) {
    mat_temp <- matrix(c(1, k * min_len, k * min_len,
                         1 + step, step + k * min_len, k * min_len),
                       ncol = 3, byrow = TRUE)
    mat_base <- rbind(mat_base, mat_temp)
    if (nbr > 1) {
      for (jj in seq_len(nbr - 1)) {
        mat_temp <- matrix(c(1 + jj * min_len, (jj + k) * min_len,
                             k * min_len, 1 + step + jj * min_len,
                             step + (jj + k) * min_len, k * min_len),
                           ncol = 3, byrow = TRUE)
        mat_base <- rbind(mat_base, mat_temp)
      }
    }
    mat_base <- mat_base[mat_base[, 2] <= t, ]
  }

  gset <- data.frame("u" = (mat_base[, 1] + mat_base[, 2]) / (2 * t),
                     "h" = mat_base[, 3] / (2 * t))

  grid_type <- "default"

  u_grid <- unique(gset[, "u"])
  h_grid <- unique(gset[, "h"])
  gset_full <- expand.grid(u = u_grid, h = h_grid)
  n         <- dim(gset)[1]
  pos_full  <- rep(TRUE, n)

  bws  <- unique(gset[, 2])
  lengths_u <- rep(NA, length(bws))
  for (i in seq_len(length(bws)))
    lengths_u[i] <- sum(gset[, 2] == bws[i])

  return(list(gset = gset, bws = h_grid, lens = lengths_u, gtype = grid_type,
              gset_full = gset_full, pos_full = pos_full))
}
