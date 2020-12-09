#' Compute Davoli immune signature
#'
#' \code{compute_davoliIS} computes Davoli immune signature as the arithmetic mean of cytotoxic
#' immune infiltrate signature genes, after rank normalization (Davoli et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Davoli immune signature
#'
#' @export
#'
#' @€xamples
#' # TODOTODO
compute.Davoli_IS <- function(RNA.tpm){

  # Literature genes
  Davoli_IS.read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")
  match_Davoli_IS.genes <- match(Davoli_IS.read, rownames(RNA.tpm))

  if (anyNA(match_Davoli_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(Davoli_IS.read[!Davoli_IS.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Davoli_IS.genes <- stats::na.omit(match_Davoli_IS.genes)
  }

  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[match_Davoli_IS.genes, ]

  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)

  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)

  message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}
