#' Compute cytolytic activity (CYT) score
#'
#' Computes CYT score as the geometric mean of its signature genes
#'
#' @references Rooney, M.S., Shukla, S.A., Wu, C.J., Getz, G., and Hacohen, N. (2015). Molecular and genetic properties
#' of tumors associated with local immune cytolytic activity. Cell 160, 48–61. https://doi.org/10.1016/j.cell.2014.12.033.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute cytolytic activity (Rooney et al, Cell, 2015)
#' CYT <- compute_CYT(RNA_tpm = gene_tpm)
compute_CYT <- function(RNA_tpm,
                        verbose = TRUE) {

  # Literature signature
  sig_read <- c("GZMA", "PRF1")
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Subset RNA_tpm
  sub_RNA_tpm <- RNA_tpm[match_sig_read, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(sub_RNA_tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  if (verbose) message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}
