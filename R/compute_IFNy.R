#' Compute IFNy signature score
#'
#' `compute_IFNy` computes IFNy signature score as the arithmetic mean of genes
#' included in the IFN-γ signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=IFNy signature score
#'
#' @export
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Compute IFNy signature (Ayers et al., JCI, 2017)
#' IFNy <- compute_IFNy(RNA_tpm = gene_tpm)
compute_IFNy <- function(RNA_tpm,
                         verbose = TRUE) {

  # Literature signature
  sig_read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA")
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", sig_read[!sig_read %in% rownames(RNA_tpm)], "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]

  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2_RNA_tpm, 2, mean)

  if (verbose) message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}
