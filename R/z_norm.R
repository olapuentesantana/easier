#' Calculation of a matrix z-score from scratch.
#'
#' \code{standarization} implements the z-score normalization.
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
standarization <- function(X, mean, sd){

  X.scale <- matrix(0, nrow(X), ncol(X), dimnames = list(rownames(X),colnames(X)))

  if (missing(mean) & missing(sd)) {
     mean.X <- colMeans(X, na.rm = T)
     sd.X <- colSds(as.matrix(X), na.rm = T)
     X.scale <- sweep(X, 2, mean.X, FUN = "-")
     X.scale <- sweep(X.scale, 2, sd.X, FUN = "/")
  } else {
    mean <- mean[na.omit(match(colnames(X), names(mean)))]
    sd <- sd[na.omit(match(colnames(X), names(sd)))]

    X <- X[, na.omit(match(names(sd), colnames(X)))]
    X.scale <- sweep(X, 2, mean, FUN = "-")
    X.scale <- sweep(X.scale, 2, sd, FUN = "/")
  }
  return(as.matrix(X.scale))
}

