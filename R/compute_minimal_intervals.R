#' Computes the set of minimal intervals as described in Duembgen (2002)
#'
#' @description Given a set of intervals, this function computes
#'              the corresponding subset of minimal intervals which are defined
#'              as follows. For a given set of intervals \eqn{\mathcal{K}},
#'              all intervals \eqn{\mathcal{I}_k \in \mathcal{K}}
#'              such that  \eqn{\mathcal{K}} does not contain a proper subset of
#'              \eqn{\mathcal{I}_k} are called minimal.
#'
#'              This function is needed for illustrative purposes.
#'              The set of all the intervals where our test rejects the null
#'              hypothesis may be quite large, hence, we would like to focus
#'              our attention on the smaller subset, for which we are still
#'              able to make simultaneous confidence intervals. This subset
#'              is the subset of minimal intervals, and it helps us to
#'              to precisely locate the intervals of further interest.
#'
#'              More details can be found in Duembgen (2002) and
#'              Khismatullina and Vogt (2019, 2020)
#' @export
#'
#' @param dataset    Set of the intervals.
#'                   It needs to contain the following columns:
#'                   "startpoint" - left end of the interval;
#'                   "endpoint"   - right end of the interval.
#' @return           Subset of minimal intervals
#'
#' @examples
#' startpoint   <- c(0, 0.5, 1)
#' endpoint     <- c(2, 2, 2)
#' dataset      <- data.frame(startpoint, endpoint)
#' minimal_ints <- compute_minimal_intervals(dataset)

compute_minimal_intervals <- function(dataset) {
  contains <- NULL
  set_cardinality <- nrow(dataset)
  dataset$endpoint <- round(dataset$endpoint, 10)
  dataset$startpoint <- round(dataset$startpoint, 10)
  if (set_cardinality > 1) {
    #Ordering the dataset such that we don't need to check previous intervals
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint), ]
    #We restore the indices after ordering
    rownames(dataset) <- seq_len(nrow(dataset))
    dataset$contains  <- numeric(set_cardinality)
    for (i in seq_len(set_cardinality - 1)) {
      for (j in (i + 1):set_cardinality) {
        if ((dataset$startpoint[i] <= dataset$startpoint[j]) &
            (dataset$endpoint[i] >= dataset$endpoint[j])) {
          #We mark all the intervals that contain at least 1 another interval.
          #We will delete them from the set afterwards.
          dataset$contains[i] <- 1
          break
        }
      }
    }
    #Now we subset everything that is not marked
    p_t_set <- subset(dataset, contains == 0, select = -contains)
  } else {
    p_t_set <- dataset
  }
  return(p_t_set)
}
