#' Compute immune cell fractions
#'
#' `compute_cell_fractions` estimates cell fractions from bulk RNA-seq data.
#'
#' Compute cell fractions from transcriptomics data. This function computes cell
#' fractions from tpm RNA-seq data using quanTIseq method
#'
#' @export
#'
#' @param RNA_tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return Cell fractions matrix: matrix of normalized enrichment scores with
#' rows=samples and columns=TFs
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(cell_fractions)
compute_cell_fractions <- function(RNA_tpm,
                                   verbose = TRUE
                                   # TODOTODO; do we need an ellipsis here?
) {

  # HGNC symbols are required
  if (any(grep("ENSG00000", rownames(RNA_tpm)))) stop("hgnc gene symbols are required", call. = FALSE)

  # Cell fractions: run deconvolute
  cell_fractions <- immunedeconv::deconvolute(gene_expression = RNA_tpm, method = "quantiseq", tumor = TRUE)

  # Samples as rows, immune cells as columns
  old_cellnames <- cell_fractions$cell_type
  new_cellnames <- c("B", "M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")

  cell_fractions <- t(cell_fractions[, -1])
  colnames(cell_fractions) <- new_cellnames
  cell_fractions[,"CD4 T"] <- cell_fractions[,"CD4 T"] +  cell_fractions[,"Treg"]

  if (verbose) message("Cell fractions computed \n")
  return(cell_fractions)
}
