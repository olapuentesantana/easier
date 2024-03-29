#' Compute Davoli immune signature (Davoli_IS) score
#'
#' Calculates Davoli_IS score as the average of the expression
#' of its signature genes after applying rank normalization, as
#' defined in Davoli et al., Science, 2017.
#'
#' @references Davoli, T., Uno, H., Wooten, E.C., and Elledge,
#' S.J. (2017). Tumor aneuploidy correlates with markers of
#' immune evasion and with reduced response to immunotherapy.
#' Science 355. https://doi.org/10.1126/science.aaf8399.
#'
#' @param matches numeric vector indicating the index of signature
#' genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols
#' in rows and samples in columns.
#'
#' @return A numeric matrix with samples in rows and Davoli_IS score
#' in a column.
#'
compute_Davoli_IS <- function(matches, RNA_tpm) {

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[matches, ]

  # Calculate rank position for each gene across samples
  ranks_sub_log2_RNA_tpm <- apply(sub_log2_RNA_tpm, 1, rank)

  # Get normalized rank by divided
  ranks_sub_log2_RNA_tpm_norm <- (ranks_sub_log2_RNA_tpm - 1) / (nrow(ranks_sub_log2_RNA_tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2_RNA_tpm_norm, 1, mean)

  return(data.frame(Davoli_IS = score, check.names = FALSE))
}
