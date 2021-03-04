#' Compute Roh immune score
#'
#' `compute_Roh_IS` computes Roh immune score as the geometric-mean of immune
#' score genes (Roh et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=Roh immune score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Roh_IS <- compute_Roh_IS(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(Roh_IS)
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

  if (verbose) message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}
