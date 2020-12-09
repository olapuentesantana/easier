#' compute_CCpair_score
#'
#' \code{compute_CCpair_score} computes the score for each cell - cell pair for all the patients.
#'
#' @export
#'
#' @param celltype1 string character with first cell type involved in the interaction
#' @param celltype2 string character with second cell type involved in the interaction
#' @param intracell.network matrix with data on cell types interaction network
#' @param lrpairs.binary binary vector displaying LR pairs with non-zero frequency
#' @param lr.frequency numeric vector with LR pairs frequency across the whole TCGA database
#' @param compute.log boolean variable in order to take the log of the weighted score
#' @param cancertype string character
#'
#' @return numeric vector with weighted scores
#'
compute_CCpair_score <- function(celltype1, celltype2, intercell.network, lrpairs.binary, lr.frequency, compute.log=T) {

  # consider the LR interactions between the two cell types
  CC.network <- intercell.network[intersect(which(intercell.network$cell1==celltype1), which(intercell.network$cell2==celltype2)),]
  CC.LRpairs <- paste(CC.network$ligands, CC.network$receptors, sep = "_")

  # extract the corresponding data for all patients
  ix <- match(CC.LRpairs, colnames(lrpairs.binary))
  CC.LR.data <- lrpairs.binary[,ix[!is.na(ix)]]

  # and the LR frequecies
  CC.LR.frequency <- lr.frequency[colnames(CC.LR.data)]

  # multiply each row of the matrix (i.e. each patient data) for the vector with the frequencies
  CC.LR.data.weighted <- t(t(CC.LR.data) * 1/CC.LR.frequency)

  # compute the cell cell interaction score as the sum of the LR weighted pairs
  CC.score <- apply(CC.LR.data.weighted, 1, sum)

  # if we use the weighted score taking the log might be better
  if (compute.log==T){
    CC.score <- log2(CC.score + 1)
  }
}
