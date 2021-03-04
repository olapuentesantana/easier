#' compute_CC_pairs_grouped
#'
#' Compute CC pairs grouped from transcriptomics data.
#'
#' `compute_CC_pairs_grouped` computes CC pairs (considering cell groups instead
#' of individual cell types) from tpm, using the null model for CC interaction
#' computed on TCGA data.
#'
#' @export
#'
#' @param lrpairs Ligand-receptor pairs weights matrix
#' @param cancertype string character
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return Cell-Cell interaction scores matrix with rows=samples and columns=cell-cell interactions
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#' lrpairs_weights <- compute_LR_pairs(
#'   RNA_tpm = Riaz_data$tpm_RNAseq,
#'   remove_genes_ICB_proxies = FALSE,
#'   cancertype = "pancan")
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs_grouped(
#'   lrpairs = lrpairs_weights,
#'   cancertype = "pancan")
#' head(ccpair_scores)
compute_CC_pairs_grouped <- function(lrpairs,
                                     cancertype = "pancan",
                                     verbose = TRUE) {

  # remove ligand receptor pairs that are always NA
  na_lrpairs <- apply(lrpairs, 2, function(x) {
    all(is.na(x))
  })
  lrpairs <- lrpairs[, na_lrpairs == FALSE]

  # binarize the data: set a threshold to 10 TPM, only pairs where both ligand and receptor have
  # TPM > 10 are kept
  lrpairs_binary <- ifelse(lrpairs > log2(10 + 1), 1, 0)
  # sum(LR.pairs.binary)/(ncol(LR.pairs.binary)*nrow(LR.pairs.binary))*100 #percentage of

  # keep only the LR.pairs for which I have (non-zero) frequencies in the TCGA
  lrpairs_binary <- lrpairs_binary[, colnames(lrpairs_binary) %in% names(lr_frequency)]

  # cancer type specific network
  intercell_network <- intercell_network_cancer_spec[[cancertype]]

  # compute the CC score for each patient
  celltypes <- unique(c(as.character(intercell_network$cell1), as.character(intercell_network$cell2)))

  CC_pairs_score <- do.call(cbind, lapply(celltypes, function(celltype1) {
    do.call(cbind, lapply(celltypes, function(celltype2) {
      compute_CCpair_score(celltype1, celltype2, intercell_network,
        lrpairs_binary, lr_frequency,
        compute_log = TRUE
      )
    }))
  }))

  metadata_CC_pairs <- do.call(rbind, lapply(celltypes, function(celltype1) {
    do.call(rbind, lapply(celltypes, function(celltype2) {
      data.frame(
        CCpair = gsub(" ", "", paste(celltype1, celltype2, sep = "_")),
        celltype1 = celltype1, celltype2 = celltype2
      )
    }))
  }))

  colnames(CC_pairs_score) <- metadata_CC_pairs$CCpair

  if (verbose) message("CC pairs computed \n")
  return(as.data.frame(CC_pairs_score))
}
