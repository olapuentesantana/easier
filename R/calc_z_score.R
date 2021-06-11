#' Z-score normalization
#'
#' Performs z-score normalization on a numeric matrix per column.
#' The user can specify mean and sd values to be used to calculate z-score values,
#' otherwise the mean and sd is calculated based on the input matrix.
#'
#' @importFrom matrixStats colSds
#'
#' @export
#'
#' @param X numeric matrix.
#' @param mean numeric vector with mean values.
#' @param sd numeric vector with sd values.
#'
#' @return Numeric matrix with values as z-scores.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # apply z-score normalization
#' tpm_zscore <- calc_z_score(t(gene_tpm))
calc_z_score <- function(X,
                         mean,
                         sd) {
  X_scale <- matrix(0, nrow(X), ncol(X), dimnames = list(rownames(X), colnames(X)))

  if (missing(mean) & missing(sd)) {
    mean_X <- colMeans(X, na.rm = TRUE)
    sd_X <- matrixStats::colSds(as.matrix(X), na.rm = TRUE)
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
