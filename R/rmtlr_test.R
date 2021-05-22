#' Function to apply regularized multi-task linear regression model
#'
#' `rmtlr_test` computes the matrix multiplication using the input features values and the parameters values derived from model optimization.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param x_test A numeric matrix containing features values (rows = samples; columns = features).
#' @param coef_matrix A numeric matrix containing the parameters values derived from model training
#' (rows = features; columns = tasks).
#'
#' @return A numeric matrix of predicted values (rows = samples; columns = tasks).
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
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = gene_tpm)
#'
#' # Parameters valus: rows = features, columns = tasks
#' estimated_parameters <- matrix(rnorm(n= (ncol(cell_fractions)+1) * 10),
#' nrow = ncol(cell_fractions) + 1, ncol = 10)
#' rownames(estimated_parameters) <- c("(Intercept)", colnames(cell_fractions))
#' colnames(estimated_parameters) <- c("CYT", "Ock_IS", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
#'
#' # Compute predictions using parameters values
#' pred_test <- rmtlr_test(x_test = cell_fractions,
#'  coef_matrix = estimated_parameters)
rmtlr_test <- function(x_test,
                       coef_matrix) {

  # Keep intercept
  intercept <- as.matrix(coef_matrix)[1, ]

  # Remove intercept from coef matrix
  coef <- as.matrix(coef_matrix)[-1, , drop = FALSE]

  # match features properly
  pos <- stats::na.omit(match(colnames(x_test), rownames(coef)))
  coef <- coef[pos, , drop = FALSE]

  if (length(coef) > 1) {
    slope <- coef
    fit_pred <- t(matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test))
    + t(slope) %*% t(as.matrix(x_test[, rownames(slope)])))
  } else {
    slope <- 0
    fit_pred <- matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test[[1]]))
  }
  return(fit_pred)
}
