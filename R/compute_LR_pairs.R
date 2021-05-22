#' Compute ligand-receptor pair weights from gene expression
#'
#' This function quantifies ligand-receptor interactions in the tumor microenvironment from
#' gene expression in TPM from bulk RNA-seq data (Lapuente-Santana et al., bioRxiv, 2021), using prior knowledge
#' coming from ligand-receptor pair annotations from the database of (Ramilowski et al., Nat Commun, 2015).
#' Each ligand-receptor weight is defined as the minimum of the log2(TPM+1) expression of the ligand and the receptor.
#'
#' @importFrom stats na.exclude
#' @importFrom utils head tail
#'
#' @param RNA_tpm A data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param remove_genes_ICB_proxies a logical value indicating whether to remove signature genes involved
#' in the derivation of hallmarks of immune response.
#' @param cancer_type A string detailing the cancer type whose ligand-receptor pairs network will be used.
#' A pan-cancer network is selected by default, whose network represents the union of all
#' ligand-receptor pairs present across the 18 cancer types studied in (Lapuente-Santana et al., bioRxiv, 2021).
#' @param verbose A logical value indicating whether to display messages about the number of ligand-receptor genes found in the gene expression data provided.
#'
#' @return A matrix of weights with samples in rows and ligand-receptor pairs in columns.
#'
#' @export
#'
#' @examples
#' # Example: Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'  install.packages("BiocManager")
#'
#' BiocManager::install(c("biomaRt",
#'  "circlize",
#'  "ComplexHeatmap",
#'  "corrplot",
#'  "DESeq2",
#'  "dplyr",
#'  "DT",
#'  "edgeR",
#'  "ggplot2",
#'  "limma",
#'  "lsmeans",
#'  "reshape2",
#'  "spatstat",
#'  "survival",
#'  "plyr"))
#'
#' install.packages("Downloads/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)
#' library(IMvigor210CoreBiologies)
#'
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(
#'   RNA_tpm = gene_tpm,
#'   remove_genes_ICB_proxies = FALSE,
#'   cancer_type = "pancan")
#' lrpair_weights[1:5, 1:5]
compute_LR_pairs <- function(RNA_tpm,
                             remove_genes_ICB_proxies = FALSE,
                             cancer_type = "pancan",
                             verbose = TRUE) {

  # Gene expression data (log2 transformed)
  gene_expr <- log2(RNA_tpm + 1)
  genes <- rownames(gene_expr)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE)

  # Genes to remove according to all ICB proxy's
  if (remove_genes_ICB_proxies) {
    if (verbose) message("Removing signatures genes for hallmarks of immune response \n")
    idy <- stats::na.exclude(match(cor_genes_to_remove, rownames(gene_expr)))
    gene_expr <- gene_expr[-idy, ]
  }

  gene_expr <- as.data.frame(gene_expr)

  # Cancer-specific LR pairs network
  intercell_network <- intercell_network_cancer_spec[[cancer_type]]
  LR_pairs <- unique(paste0(intercell_network$ligands, "_", intercell_network$receptors))

  # check what is the percentage of genes we have in our data
  all_lrpairs_genes <- unique(c(intercell_network$ligands, intercell_network$receptors))
  genes_kept <- intersect(rownames(gene_expr), all_lrpairs_genes)
  genes_left <- setdiff(all_lrpairs_genes, rownames(gene_expr))

  # check what is the percentage of regulated transcripts that we have in our data
  message("LR signature genes found in data set: ", length(genes_kept), "/", length(all_lrpairs_genes), " (", round(length(genes_kept)/length(all_lrpairs_genes), 3) * 100,"%)")

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

  if (verbose) message("LR pairs computed \n")
  return(as.data.frame(LR.pairs.computed))
}
