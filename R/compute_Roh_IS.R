#' Compute Roh immune score (Roh_IS)
#'
#' This function computes Roh_IS score as the geometric-mean of its signature genes.
#'
#' @references Roh, W., Chen, P.-L., Reuben, A., Spencer, C.N., Prieto, P.A., Miller, J.P., Gopalakrishnan, V.,
#' Wang, F., Cooper, Z.A., Reddy, S.M., et al. (2017). Integrated molecular analysis of tumor biopsies on sequential
#' CTLA-4 and PD-1 blockade reveals markers of response and resistance. Sci. Transl. Med. 9.
#' https://doi.org/10.1126/scitranslmed.aah3560.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and Roh_IS score in a column.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute Roh immune score (Roh et al., Sci. Transl. Med., 2017)
#' Roh_IS <- compute_Roh_IS(RNA_tpm = gene_tpm)
compute_Roh_IS <- function(RNA_tpm,
                           verbose = TRUE) {

  # Literature signature
  sig_read <- c(
    "GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
    "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
    "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
    "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4",
    "ICAM5", "VCAM1"
  )

  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Subset RNA_tpm
  sub_RNA_tpm <- RNA_tpm[match_sig_read, ]

  # Pseudocount of 0.01 for all genes
  sub_RNA_tpm <- sub_RNA_tpm + 0.01

  # Pseudocount of 1 for genes with 0 expr
  if (any(sub_RNA_tpm == 0)) sub_RNA_tpm[sub_RNA_tpm == 0] <- sub_RNA_tpm[sub_RNA_tpm == 0] + 1

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_RNA_tpm, 2, function(X) exp(mean(log(X))))

  if (verbose) message("Roh_IS score computed")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}
