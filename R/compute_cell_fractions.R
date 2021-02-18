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
#'
#' @return Cell fractions matrix: matrix of normalized enrichment scores with
#' rows=samples and columns=TFs
#'
#' @examples
#' # TODOTODO
compute_cell_fractions <- function(RNA_tpm
                                   # TODOTODO; do we need an ellipsis here?
                                   ){

  # ****************
  # packages

  # TODOTODO: we should handle this part outside the function call, i.e. in the dependencies - might require we go fully fledged with immunedeconv
  # if(!("BiocManager" %in% installed.packages()[,"Package"])) install.packages("BiocManager", quiet = TRUE)
  # list.of.packages <- c("remotes")
  # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  # if(length(new.packages)) BiocManager::install(new.packages, ask = FALSE)
  #
  # suppressMessages(remotes::install_github("icbi-lab/immunedeconv"))

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", rownames(RNA_tpm)))) stop("hgnc gene symbols are required", call. = FALSE))

  # Cell fractions: run deconvolute
  cell_fractions <- immunedeconv::deconvolute(gene_expression = RNA_tpm, method = "quantiseq", tumor = TRUE)

  # Samples as rows, TFs as columns
  old_cellnames <- cell_fractions$cell_type
  new_cellnames <- c("B","M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")

  cell_fractions <- t(cell_fractions[,-1])
  colnames(cell_fractions) <- new_cellnames

  message("Cell fractions computed \n")

  return(cell_fractions)
}
