#' Compute chemokine score
#'
#' \code{compute_chemokine} computes chemoine score as the PC1 score that results from applying PCA
#' to z-score expression of 12 chemokine genes (Messina et al., 2012).
#'
#' @importFrom stats na.omit prcomp
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
#'
#' @€xamples
#' # TODOTODO
compute.chemokines <- function(RNA.tpm){

  # Literature genes
  chemokines.read <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                      "CXCL9", "CXCL10", "CXCL11", "CXCL13")
  match_chemokines.genes <- match(chemokines.read, rownames(RNA.tpm))

  if (anyNA(match_chemokines.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(chemokines.read[!chemokines.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_chemokines.genes <- stats::na.omit(match_chemokines.genes)
  }

  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset gene_expr
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_chemokines.genes, ]

  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine.pca <- stats::prcomp(t(sub_log2.RNA.tpm), center = T, scale = T)
  score <- chemokine.pca$x[, 1]

  message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}
