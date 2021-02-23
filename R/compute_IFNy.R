#' Compute IFNy signature score
#'
#' `compute_IFNy` computes IFNy signature score as the arithmetic mean of genes
#' included in the IFN-γ signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IFNy signature score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' IFNy <- compute_IFNy(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(IFNy)
compute_IFNy <- function(RNA_tpm) {

  # Literature genes
  IFNy.read <- c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA")
  match_IFNy.genes <- match(IFNy.read, rownames(RNA_tpm))

  if (anyNA(match_IFNy.genes)) {
    warning(paste0("differenty named or missing signature genes : \n", IFNy.read[!IFNy.read %in% rownames(RNA_tpm)]))
    match_IFNy.genes <- stats::na.omit(match_IFNy.genes)
  }

  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2.RNA_tpm <- log2.RNA_tpm[match_IFNy.genes, ]

  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2.RNA_tpm, 2, mean)

  message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}
