#' Compute Expanded Immune signature (Ayers_expIS) score
#'
#' This function calculates Ayers_expIS score as the average expression of its signature genes.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with rows=samples and columns=Expanded Immune signature score.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute expanded immune signature (Ayers et al., JCI, 2017)
#' Ayers_expIS <- compute_Ayers_expIS(gene_tpm)
compute_Ayers_expIS <- function(RNA_tpm,
                                verbose = TRUE) {

  # Literature signature
  sig_read <- c(
    "GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1", "CD3D", "CD3E",
    "CD2", "IL2RG", "NKG7", "HLA-E", "CIITA", "HLA-DRA", "LAG3", "IDO1", "TAGAP"
  )
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differenty named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2_RNA_tpm
  sub_log2_RNA_tpm <- log2_RNA_tpm[match_sig_read, ]

  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2_RNA_tpm, 2, mean)

  if (verbose) message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}
