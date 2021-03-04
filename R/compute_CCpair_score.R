#' compute_CCpair_score
#'
#' `compute_CCpair_score` computes the score for each cell-cell pair for all the patients.
#'
#' @export
#'
#' @param celltype1 string character with first cell type involved in the interaction
#' @param celltype2 string character with second cell type involved in the interaction
#' @param intercell_network matrix with data on cell types interaction network
#' @param lrpairs_binary binary vector displaying LR pairs with non-zero frequency
#' @param lr_frequency numeric vector with LR pairs frequency across the whole TCGA database
#' @param compute_log boolean variable in order to take the log of the weighted score
#' @param cancertype string character
#'
#' @return numeric vector with weighted scores
#'
#' @examples
#' # TODOTODO
compute_CCpair_score <- function(celltype1,
                                 celltype2,
                                 intercell_network,
                                 lrpairs_binary,
                                 lr_frequency,
                                 compute_log = TRUE,
                                 cancertype) {

  # consider the LR interactions between the two cell types
  CC_network <- intercell_network[intersect(which(intercell_network$cell1 == celltype1), which(intercell_network$cell2 == celltype2)), ]
  CC_LRpairs <- paste(CC_network$ligands, CC_network$receptors, sep = "_")

  # extract the corresponding data for all patients
  ix <- match(CC_LRpairs, colnames(lrpairs_binary))
  CC_LR_data <- lrpairs_binary[, ix[!is.na(ix)]]

  # and the LR frequecies
  CC_lr_frequency <- lr_frequency[colnames(CC_LR_data)]

  # multiply each row of the matrix (i.e. each patient data) for the vector with the frequencies
  CC_LR_data_weighted <- t(t(CC_LR_data) * 1 / CC_lr_frequency)

  # compute the cell cell interaction score as the sum of the LR weighted pairs
  CC_score <- apply(CC_LR_data_weighted, 1, sum)

  # if we use the weighted score taking the log might be better
  if (compute_log == TRUE) {
    CC_score <- log2(CC_score + 1)
  }
}
