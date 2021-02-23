#' Compute cytolytic activity score
#'
#' `compute_CYT` computes cytolytic activity score as the geometric mean of immune
#' cytolytic genes (Rooney et al., 2015).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' CYT <- compute_CYT(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(CYT)
compute_CYT <- function(RNA_tpm) {

  # Literature genes
  CYT.read <- c("GZMA", "PRF1")
  match_CYT.genes <- match(CYT.read, rownames(RNA_tpm))

  if (anyNA(match_CYT.genes)) {
    warning(paste0("differenty named or missing signature genes : \n", paste(CYT.read[!CYT.read %in% rownames(RNA_tpm)], collapse = "\n")))
    match_CYT.genes <- stats::na.omit(match_CYT.genes)
  }

  # Subset RNA_tpm
  subset_RNA_tpm <- RNA_tpm[match_CYT.genes, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA_tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}
