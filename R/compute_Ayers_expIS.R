#' Compute Expanded Immune signature
#'
#' `compute_Ayers_expIS` computes Expanded Immune signature score as the arithmetic
#' mean of genes included in the Expanded Immune signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Expanded Immune signature score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Ayers_expIS <- compute_Ayers_expIS(RNA_tpm= Riaz_data$tpm_RNAseq)
compute_Ayers_expIS <- function(RNA_tpm) {

  # Literature genes
  Ayers_expIS.read <- c(
    "GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1", "CD3D", "CD3E",
    "CD2", "IL2RG", "NKG7", "HLA-E", "CIITA", "HLA-DRA", "LAG3", "IDO1", "TAGAP"
  )
  match_Ayers_expIS.genes <- match(Ayers_expIS.read, rownames(RNA_tpm))

  if (anyNA(match_Ayers_expIS.genes)) {
    warning(c("differenty named or missing signature genes : \n", paste(Ayers_expIS.read[!Ayers_expIS.read %in% rownames(RNA_tpm)], collapse = "\n")))
    match_Ayers_expIS.genes <- stats::na.omit(match_Ayers_expIS.genes)
  }

  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2.RNA_tpm <- log2.RNA_tpm[match_Ayers_expIS.genes, ]

  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2.RNA_tpm, 2, mean)

  message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}
