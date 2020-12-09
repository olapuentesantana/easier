#' Compute the expression of the immune checkpoints genes
#'
#' \code{computation_ICB_genes} computes the scores for the immune checkpoint genes.
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#'
#' @return List with the expression of the immune checkpoint genes
#'
#' @€xamples
#' # TODOTODO
compute_ICB_genes <- function(RNA.tpm){

  # Extract position genes for GZMA and PRF1
  tmp <- match(c("CD274","CTLA4","PDCD1"), rownames(RNA.tpm))

  # PDL-1 calculation
  PDL1_expr = RNA.tpm[tmp[1],] ; rownames(PDL1_expr) <- "PDL1"

  # CTLA-4 calculation
  CTLA4_expr = RNA.tpm[tmp[2],] ; rownames(CTLA4_expr) <- "CTLA4"

  # PD-1 calculation
  PD1_expr = RNA.tpm[tmp[3],] ; rownames(PD1_expr) <- "PD1"

  ICB_genes_expr <- list(PDL1 = PDL1_expr , CTLA4 = CTLA4_expr , PD1 =PD1_expr)
  message("ICB genes expression computed")
  return(ICB_genes_expr)
}
