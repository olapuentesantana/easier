#' Compute ligand-receptor pairs
#'
#' \code{compute_LR_pairs} obtain ligand-receptor pairs weights from tpm RNAseq data.
#'
#' @importFrom stats na.exclude
#' @importFrom utils head tail
#'
#' @param RNA.tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove.genes.ICB_proxies boolean variable to reomove all those genes involved in the computation of ICB proxy's of response
#' @param compute.cytokines.pairs boolean variable to compute cytokine pairs as well
#' @param cancertype string character
#'
#' @return A list with the following elements:
#'         \describe{
#'               \item{LRpairs}{Ligand-leceptor pairs weights matrix in log2(tpm + 1) with rows=samples and columns = L-R pairs}
#'               \item{CYTOKINEpairs}{Cytokine-Cytokine pairs weights matrix in log2(tpm + 1) with rows=samples and columns = CYTOKINE pairs}
#'         }
#' @export
#'
#' @examples
#' # TODOTODO
compute_LR_pairs <- function(RNA.tpm,
                             remove.genes.ICB_proxies=F,
                             compute.cytokines.pairs=F,
                             cancertype){

  # Gene expression data (log2 transformed)
  gene_expr <- log2(RNA.tpm + 1)
  genes <- rownames(gene_expr)

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE))

  # Genes to remove according to all ICB proxy's
  if (remove.genes.ICB_proxies) {
    message("Removing signatures genes for proxy's of ICB response  \n")
    idy <- stats::na.exclude(match(all_genes_to_remove, rownames(gene_expr)))
    gene_expr <- gene_expr[-idy,]
  }

  gene_expr <- as.data.frame(gene_expr)

  # Cancer-specific LR pairs network
  intercell.network <- intercell.network.cancer.spec[[cancertype]]
  LR_pairs <- unique(paste0(intercell.network$ligands, "_", intercell.network$receptors))

  # Compute L-R pairs
  LR.pairs.computed <- do.call(rbind, lapply(1:length(LR_pairs), function(x){

    ligand <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), head, 1)
    receptor <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), tail, 1)

    pos_lr <- match(c(ligand, receptor), rownames(gene_expr))
    # When a ligand or receptor is not found, NA value should be returned.
    by_patient <- t(as.data.frame(apply(gene_expr[pos_lr, ], 2, min)))
    rownames(by_patient) <- LR_pairs[x]
    return(by_patient)
  }))
  LR.pairs.computed <- t(LR.pairs.computed)

  # Compute cytokine pairs
  if (compute.cytokines.pairs) {
    idy <- stats::na.exclude(match(CYTOKINE.pairs_subnetwork, colnames(LR.pairs.computed)))
    CYTOKINE.pairs.computed <- LR.pairs.computed[,idy]

    # Apply grouping to LRpairs data
    for (X in 1:length(grouping_lrpairs_info)){

      keep <- unique(grouping_lrpairs_info[[X]]$main)
      remove <- unique(grouping_lrpairs_info[[X]]$involved_pairs)
      combo_name <- unique(grouping_lrpairs_info[[X]]$combo_name)

      pos_remove <- match(remove, colnames(LR.pairs.computed))
      pos_keep <- match(keep, colnames(LR.pairs.computed))

      colnames(LR.pairs.computed)[pos_keep] <- combo_name
      LR.pairs.computed <- LR.pairs.computed[, -pos_remove]

      pos_remove <- match(remove, colnames(CYTOKINE.pairs.computed))
      pos_keep <- match(keep, colnames(CYTOKINE.pairs.computed))

      if(all(is.na(pos_remove) == F) & all(is.na(pos_keep) == F)){
        colnames(CYTOKINE.pairs.computed)[pos_keep] <- combo_name
        CYTOKINE.pairs.computed <- CYTOKINE.pairs.computed[, -pos_remove]
      }
    }
    LR.data <- list(LRpairs = as.data.frame(LR.pairs.computed), CYTOKINEpairs = as.data.frame(CYTOKINE.pairs.computed))
  }else{

    # Apply grouping to LRpairs data
    for (X in 1:length(grouping_lrpairs_info)){

      keep <- unique(grouping_lrpairs_info[[X]]$main)
      remove <- unique(grouping_lrpairs_info[[X]]$involved_pairs)
      combo_name <- unique(grouping_lrpairs_info[[X]]$combo_name)

      pos_remove <- match(remove, colnames(LR.pairs.computed))
      pos_keep <- match(keep, colnames(LR.pairs.computed))

      colnames(LR.pairs.computed)[pos_keep] <- combo_name
      LR.pairs.computed <- LR.pairs.computed[, -pos_remove]
    }
    LR.data <- list(LRpairs = as.data.frame(LR.pairs.computed))
  }

  message("L-R pairs computed \n")
  message("Cytokine pairs computed \n")
  return(LR.data)
}

