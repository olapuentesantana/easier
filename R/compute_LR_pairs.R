#' Compute ligand-receptor pairs
#'
#' `compute_LR_pairs` obtain ligand-receptor pairs weights from tpm RNA-seq data.
#'
#' @importFrom stats na.exclude
#' @importFrom utils head tail
#'
#' @param RNA_tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove_genes_ICB_proxies boolean variable to remove all those genes
#' involved in the computation of ICB proxy's of response
#' @param cancertype string character
#'
#' @return Ligand-receptor pairs weights matrix in log2(tpm + 1) with rows=samples and columns=L-R pairs
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' # Computation of ligand-receptor pair weights
#' lrpairs_weights <- compute_LR_pairs(
#'   RNA_tpm = Riaz_data$tpm_RNAseq,
#'   remove_genes_ICB_proxies = FALSE,
#'   cancertype = "pancan")
#' lrpairs_weights[1:5, 1:5]
compute_LR_pairs <- function(RNA_tpm,
                             remove_genes_ICB_proxies = FALSE,
                             cancertype = "pancan") {

  # Gene expression data (log2 transformed)
  gene_expr <- log2(RNA_tpm + 1)
  genes <- rownames(gene_expr)

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE))

  # Genes to remove according to all ICB proxy's
  if (remove_genes_ICB_proxies) {
    message("Removing signatures genes for proxy's of ICB response  \n")
    idy <- stats::na.exclude(match(cor_genes_to_remove, rownames(gene_expr)))
    gene_expr <- gene_expr[-idy, ]
  }

  gene_expr <- as.data.frame(gene_expr)

  # Cancer-specific LR pairs network
  intercell_network <- intercell_network_cancer_spec[[cancertype]]
  LR_pairs <- unique(paste0(intercell_network$ligands, "_", intercell_network$receptors))

  # Compute L-R pairs
  LR.pairs.computed <- do.call(rbind, lapply(1:length(LR_pairs), function(x) {
    ligand <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), head, 1)
    receptor <- sapply(strsplit(LR_pairs[x], split = "_", fixed = TRUE), tail, 1)

    pos_lr <- match(c(ligand, receptor), rownames(gene_expr))
    # When a ligand or receptor is not found, NA value should be returned.
    by_patient <- t(as.data.frame(apply(gene_expr[pos_lr, ], 2, min)))
    rownames(by_patient) <- LR_pairs[x]
    return(by_patient)
  }))
  LR.pairs.computed <- t(LR.pairs.computed)

  # Apply grouping to LRpairs data
  for (X in 1:length(grouping_lrpairs_info)) {
    keep <- unique(grouping_lrpairs_info[[X]]$main)
    remove <- unique(grouping_lrpairs_info[[X]]$involved_pairs)
    combo_name <- unique(grouping_lrpairs_info[[X]]$combo_name)

    pos_remove <- match(remove, colnames(LR.pairs.computed))
    pos_keep <- match(keep, colnames(LR.pairs.computed))

    colnames(LR.pairs.computed)[pos_keep] <- combo_name
    LR.pairs.computed <- LR.pairs.computed[, -pos_remove]
  }

  message("LR pairs computed \n")
  return(as.data.frame(LR.pairs.computed))
}
