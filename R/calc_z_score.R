#' Z-score normalization
#'
#' `calc_z_score` applies z-score normalization on a numeric matrix per column.
#' The user can specify mean and sd values to be used to calculate z-score values,
#' otherwise the mean and sd is calculated based on the input matrix.
#'
#' @importFrom matrixStats colSds
#'
#' @export
#'
#' @param X A numeric matrix.
#' @param mean A numeric vector with mean values.
#' @param sd A numeric vector with sd values.
#'
#' @return A numeric matrix with values as z-scores.
#'
#' @examples
#' # Example: Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'  install.packages("BiocManager")
#'
#' BiocManager::install(c("biomaRt",
#'  "circlize",
#'  "ComplexHeatmap",
#'  "corrplot",
#'  "DESeq2",
#'  "dplyr",
#'  "DT",
#'  "edgeR",
#'  "ggplot2",
#'  "limma",
#'  "lsmeans",
#'  "reshape2",
#'  "spatstat",
#'  "survival",
#'  "plyr"))
#'
#' install.packages("Downloads/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)
#' library(IMvigor210CoreBiologies)
#'
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' tpm_zscore <- calc_z_score(t(gene_tpm))
#' head(tpm_zscore)
calc_z_score <- function(X,
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
