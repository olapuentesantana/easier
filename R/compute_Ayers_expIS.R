#' Compute Expanded Immune signature
#'
#' `compute_Ayers_expIS` computes Expanded Immune signature score as the arithmetic
#' mean of genes included in the Expanded Immune signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=Expanded Immune signature score
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
#' Ayers_expIS <- compute_Ayers_expIS(gene_tpm)
#' head(Ayers_expIS)
compute_Ayers_expIS <- function(RNA_tpm,
                                verbose = TRUE) {

  # Literature signature
  sig_read <- c(
    "GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1", "CD3D", "CD3E",
    "CD2", "IL2RG", "NKG7", "HLA-E", "CIITA", "HLA-DRA", "LAG3", "IDO1", "TAGAP"
  )
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2_RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]

  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2_RNA_tpm, 2, mean)

  if (verbose) message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}
