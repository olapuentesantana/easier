#' Compute Davoli immune signature (Davoli_IS) score
#'
#' The function calculates Davoli_IS score as the average of the expression of its signature genes after applying rank normalization
#'
#' @references Davoli, T., Uno, H., Wooten, E.C., and Elledge, S.J. (2017). Tumor aneuploidy correlates
#' with markers of immune evasion and with reduced response to immunotherapy. Science 355.
#' https://doi.org/10.1126/science.aaf8399.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return A numeric matrix with samples in rows and Davoli_IS score in a column.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute davoli immune signature (Davoli et al., Science 2017)
#' Davoli_IS <- compute_Davoli_IS(RNA_tpm = gene_tpm)
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
