#' Calculation of a matrix z-score from scratch.
#'
#' `standardization` implements the z-score normalization.
#'
#' @importFrom matrixStats colSds
#'
#' @export
#'
#' @param X numeric matrix with data
#' @param mean numeric vector with data
#' @param sd numeric vector with data
#'
#' @return numeric matrix with scaled data
#'
#' @examples
#' # TODOTODO
standardization <- function(X,
                            mean,
                            sd) {
  X_scale <- matrix(0, nrow(X), ncol(X), dimnames = list(rownames(X), colnames(X)))

  if (missing(mean) & missing(sd)) {
    mean_X <- colMeans(X, na.rm = TRUE)
    sd_X <- colSds(as.matrix(X), na.rm = TRUE)
    X_scale <- sweep(X, 2, mean_X, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd_X, FUN = "/")
  } else {
    mean <- mean[na.omit(match(colnames(X), names(mean)))]
    sd <- sd[na.omit(match(colnames(X), names(sd)))]

    X <- X[, na.omit(match(names(sd), colnames(X)))]
    X_scale <- sweep(X, 2, mean, FUN = "-")
    X_scale <- sweep(X_scale, 2, sd, FUN = "/")
  }
  return(as.matrix(X_scale))
}
