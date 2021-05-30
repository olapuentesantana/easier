#' Compute immune cell fractions from gene expression using quanTIseq
#'
#' This function estimates cell fractions from bulk RNA-seq data (tpm)
#' using quanTIseq method from (Finotello et al., Genome Med, 2019).
#'
#' @export
#'
#' @param RNA_tpm A data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose A logical value indicating whether to display messages about the number of immune cell
#' signature genes found in the gene expression data provided.
#'
#' @return A matrix of normalized enrichment scores with samples in rows and cell types in columns.
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = gene_tpm)
compute_cell_fractions <- function(RNA_tpm,
                                   verbose = TRUE
                                   # TODOTODO; do we need an ellipsis here?
) {

  # Some checks
  if (is.null(RNA_tpm)) stop("tpm data not found")

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
