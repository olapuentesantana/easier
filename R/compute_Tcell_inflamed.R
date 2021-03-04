#' Compute T cell-inflamed signature score
#'
#' `compute_ayersTcellInfl` computes T cell-inflamed signature score by
#' taking a weighted sum of the housekeeping normalized values of the T cell-inflamed
#' signature genes
#' TODOTODOoscar, async name of function vs doc
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=T cell-inflamed signature
#' score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Tcell_inflamed <- compute_Tcell_inflamed(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(Tcell_inflamed)
compute_Tcell_inflamed <- function(RNA_tpm,
                                   verbose = TRUE) {

  # Literature signature
  sig_read <- c(
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
    "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"
  )
  house_read <- c(
    "STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2",
    "UBB", "TBP", "SDHA"
  ) # C14orf102 = NRDE2
  weights <- data.frame(
    CCL5 = 0.008346, CD27 = 0.072293, CD274 = 0.042853, CD276 = -0.0239, CD8A = 0.031021,
    CMKLR1 = 0.151253, CXCL9 = 0.074135, CXCR6 = 0.004313, `HLA-DQA1` = 0.020091,
    `HLA-DRB1` = 0.058806, `HLA-E` = 0.07175, IDO1 = 0.060679, LAG3 = 0.123895,
    NKG7 = 0.075524, PDCD1LG2 = 0.003734, PSMB10 = 0.032999, STAT1 = 0.250229,
    TIGIT = 0.084767, check.names = FALSE
  )

  # Some genes might have other name: case for "C14orf102", it's called "NRDE2", be careful
  if (any(rownames(RNA_tpm) %in% "C14orf102")) {
    warning("Gene name changed: NRDE2 is approved symbol, not C14orf102", "\n")
    rownames(RNA_tpm)[rownames(RNA_tpm) %in% "C14orf102"] <- "NRDE2"
  }

  match_house_read <- match(house_read, rownames(RNA_tpm))
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(c(match_house_read, match_sig_read))) {
    tmp <- c(sig_read, house_read)
    warning("differenty named or missing signature genes : \n", paste(tmp[!tmp %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_house_read <- stats::na.omit(match_house_read)
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2_RNA_tpm
  ## housekeeping
  house_log2_RNA_tpm <- log2_RNA_tpm[match_house_read, ]
  ## predictors
  sig_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]
  weights <- weights[, rownames(sig_log2_RNA_tpm)]

  # Housekeeping normalization
  average_house_log2_RNA_tpm <- apply(house_log2_RNA_tpm, 2, mean)
  norm_sig_log2_RNA_tpm <- sweep(sig_log2_RNA_tpm, 2, average_house_log2_RNA_tpm, FUN = "-")

  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(norm_sig_log2_RNA_tpm), colnames(as.vector(weights)))
  score <- t(norm_sig_log2_RNA_tpm[tidy, ]) %*% t(as.vector(weights))

  if (verbose) message("Tcell_inflamed score computed")
  return(data.frame(Tcell_inflamed = score, check.names = FALSE))
}
