#' Compute microsatellite instability status (MSI) score
#'
#' Computes MSI score by logical comparison of MSI-related gene pairs expression
#'
#' @references Fu, Y., Qi, L., Guo, W., Jin, L., Song, K., You, T., Zhang, S., Gu, Y., Zhao, W.,
#' and Guo, Z. (2019). A qualitative transcriptional signature for predicting microsatellite instability
#' status of right-sided Colon Cancer. BMC Genomics 20, 769. https://doi.org/10.1186/s12864-019-6129-8.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=MSI score
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute Microsatellite instability status (Fu et al., BMC Genomics, 2019 )
#' MSI <- compute_MSI(RNA_tpm = gene_tpm)
compute_MSI <- function(RNA_tpm,
                        verbose = TRUE) {

  # Literature signature: * (CCRN4L in tcga, NOCT approved symbol)
  sig_basis <- data.frame(
    Gene_1 = c(
      "HNRNPL", "MTA2", "CALR", "RASL11A", "LYG1", "STRN3", "HPSE",
      "PRPF39", "NOCT", "AMFR"
    ),
    Gene_2 = c(
      "CDC16", "VGF", "SEC22B", "CAB39L", "DHRS12", "TMEM192", "BCAS3",
      "ATF6", "GRM8", "DUSP18"
    )
  )

  sig_read <- unique(as.vector(as.matrix(sig_basis))) # 20 genes

  # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
  if (any(rownames(RNA_tpm) %in% "CCRN4L")) {
    if (verbose) warning("Gene name changed: NOCT is approved symbol, not CCRN4L", "\n")
    rownames(RNA_tpm)[rownames(RNA_tpm) %in% "CCRN4L"] <- "NOCT"
  }

  # Subset RNA_tpm
  match_F_1 <- match(as.character(sig_basis[, 1]), rownames(RNA_tpm))
  match_F_2 <- match(as.character(sig_basis[, 2]), rownames(RNA_tpm))

  if (anyNA(c(match_F_1, match_F_2))) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
  }

  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  F_pair_expr_B <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  F_matrix <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  colnames(F_matrix) <- colnames(RNA_tpm)
  remove_pairs <- vector("list", length = ncol(RNA_tpm))
  names(remove_pairs) <- colnames(RNA_tpm)
  score <- vector("numeric", length = ncol(RNA_tpm))
  names(score) <- colnames(RNA_tpm)

  # Log2 transformation:
  log2_RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2_RNA_tpm[match_F_1, ]
  F_pair_expr_B <- log2_RNA_tpm[match_F_2, ]

  if (anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }

  F_matrix <- F_pair_expr_A > F_pair_expr_B
  if (anyNA(F_matrix)) {
    score <- colSums(F_matrix, na.rm = TRUE)
    score <- (score * nrow(F_matrix)) / (nrow(F_matrix) - length(remove_pairs))
  } else {
    score <- colSums(F_matrix)
  }

  if (verbose) message("MSI score computed")
  return(data.frame(MSI = score, check.names = FALSE))
}
