#' Compute chemokine score
#'
#' `compute_chemokines` computes chemokine score as the PC1 score that results from
#' applying PCA to z-score expression of 12 chemokine genes (Messina et al., 2012).
#'
#' @importFrom stats na.omit prcomp
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
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
#' chemokines <- compute_chemokines(RNA_tpm = gene_tpm)
#' head(chemokines)
compute_chemokines <- function(RNA_tpm,
                               verbose = TRUE) {

  # Literature genes
  sig_read <- c(
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
    "CXCL9", "CXCL10", "CXCL11", "CXCL13"
  )

  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset gene_expr
  sub_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]

  # calculation: using PCA (Z-score calculated within prcomp)
  pcs <- stats::prcomp(t(sub_log2_RNA_tpm), center = TRUE, scale = TRUE)
  score <- pcs$x[, 1]

  if (verbose) message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}
