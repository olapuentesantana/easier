#' Compute T cell-inflamed signature (Tcell_inflamed) score
#'
#' This function calculates Tcell_inflamed score as a weighted sum of housekeeping normalized expression of its signature genes.
#' Weightes were available at Table S2B from Cristescu R, et al. Pan-tumor genomic biomarkers for PD-1 checkpoint
#' blockade-based immunotherapy. Science. (2018) 362:eaar3593. doi: 10.1126/science.aar3593.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return A numeric matrix with samples in rows and Tcell_inflamed score in a column.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute T-cell inflamed signature (Ayers et al., JCI, 2017)
#' Tcell_inflamed <- compute_Tcell_inflamed(RNA_tpm = gene_tpm)
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
