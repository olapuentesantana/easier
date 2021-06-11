#' Computation of tertiary lymphoid structures signature (TLS) score
#'
#' Computes TLS score as the geometric-mean of the expression of its signature genes.
#'
#' @references Cabrita, R., Lauss, M., Sanna, A., Donia, M., Skaarup Larsen, M., Mitra, S., Johansson, I., Phung, B.,
#' Harbst, K., Vallon-Christersson, J., et al. (2020). Tertiary lymphoid structures improve immunotherapy and survival
#' in melanoma. Nature 577, 561–565.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=TLS signature
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute tertiary lymphoid structures signature (Cabrita et al., Nature, 2020)
#' TLS <- compute_TLS(RNA_tpm = gene_tpm)
compute_TLS <- function(RNA_tpm,
                        verbose = TRUE) {

  # Literature signature
  sig_read <- c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS")
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Subset RNA_tpm
  sub_RNA_tpm <- RNA_tpm[match_sig_read, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_RNA_tpm, 2, function(X) exp(mean(log2(X + 1))))

  if (verbose) message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}
