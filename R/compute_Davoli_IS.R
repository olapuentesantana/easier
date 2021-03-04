#' Compute Davoli immune signature
#'
#' `compute_Davoli_IS` computes Davoli immune signature as the arithmetic mean of cytotoxic
#' immune infiltrate signature genes, after rank normalization (Davoli et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=Davoli immune signature
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Davoli_IS <- compute_Davoli_IS(RNA_tpm= Riaz_data$tpm_RNAseq)
#' head(Davoli_IS)
compute_Davoli_IS <- function(RNA_tpm,
                              verbose = TRUE) {

  # Literature signature
  sig_read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning(c("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n")))
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]

  # Calculate rank position for each gene across samples
  rank_sub_log2_RNA_tpm <- apply(sub_log2_RNA_tpm, 1, rank)

  # Get normalized rank by divided
  norm_rank_sub_log2_RNA_tpm <- (rank_sub_log2_RNA_tpm - 1) / (nrow(rank_sub_log2_RNA_tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(norm_rank_sub_log2_RNA_tpm, 1, mean)

  if (verbose) message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}
