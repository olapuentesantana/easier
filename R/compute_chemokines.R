#' Compute chemokine signature (chemokines) score
#'
#' Computes chemokines score as the PC1 score that results from applying PCA to z-score expression of its signature genes
#'
#' @references Messina, J.L., Fenstermacher, D.A., Eschrich, S., Qu, X., Berglund, A.E., Lloyd, M.C., Schell, M.J.,
#' Sondak, V.K., Weber, J.S., and Mulé, J.J. (2012). 12-Chemokine gene signature identifies lymph node-like structures
#' in melanoma: potential for patient selection for immunotherapy? Sci. Rep. 2, 765. https://doi.org/10.1038/srep00765.
#'
#' @importFrom stats na.omit prcomp
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute chemokine signature (Messina et al., Nat. Sci. Rep., 2012)
#' chemokines <- compute_chemokines(RNA_tpm = gene_tpm)
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
