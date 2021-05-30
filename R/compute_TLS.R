#' Compute tertiary lymphoid structures signature
#'
#' `compute_TLS` computes TLS signature as the geometric-mean of TLS signature
#' genes (Cabrita et al., 2020).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=TLS signature
#'
#' @export
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
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
