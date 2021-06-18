#' Compute pathway activity from gene expression using PROGENy
#'
#' Infers pathway activity from gene expression in raw counts from bulk RNA-seq data
#' using PROGENy method from (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018).
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#' getVarianceStabilizedData
#' @import progeny
#' @importFrom stats na.exclude
#'
#' @param RNA_counts data.frame containing raw counts values (with HGNC gene symbols as row names and samples identifiers as column names).
#' @param remove_genes_ICB_proxies logical value indicating whether to remove signature genes involved
#' in the derivation of hallmarks of immune response.
#' @param verbose logical value indicating whether to display messages about the number of pathway signature
#' genes found in the gene expression data provided.
#'
#' @return A matrix with samples in rows and pathways in columns.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_count <- dataset_mariathasan@counts
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018)
#' pathway_activity <- compute_pathways_scores(
#'   RNA_counts = gene_count,
#'   remove_genes_ICB_proxies = TRUE
#' )
compute_pathways_scores <- function(RNA_counts,
                                    remove_genes_ICB_proxies = TRUE,
                                    verbose = TRUE) {
  # Some checks
  if (is.null(RNA_counts)) stop("raw counts data not found")

  # Gene expression data
  raw_counts <- RNA_counts
  genes <- rownames(raw_counts)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE)

  # Remove list of genes used to build proxy's of ICB response
  if (remove_genes_ICB_proxies) {
    if (verbose) message("Removing signatures genes for proxy's of ICB response  \n")
    idy <- stats::na.exclude(match(cor_genes_to_remove, rownames(raw_counts)))
    raw_counts <- raw_counts[-idy, ]
  }

  # Integers are required for "DESeq2"
  if (is.integer(raw_counts) == FALSE) {
    raw_counts_integer <- apply(raw_counts, 2, as.integer)
    rownames(raw_counts_integer) <- rownames(raw_counts)
  } else {
    raw_counts_integer <- raw_counts
  }

  # Variance stabilizing transformation (DESeq2 package)
  # Integer count matrix, a data frame with the sample info,  design =~1 to consider all samples as part of the same group.

  # Column data:
  colData <- data.frame(id = colnames(raw_counts_integer))

  if (verbose) message("DESeq2 normalization -->\n")
  # Construction a DESeqDataSet: (Forced all to be data.frames($ operator))
  dset <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts_integer,
    colData = colData,
    design = ~1
  )

  # Variance stabilization transformation
  dset <- DESeq2::estimateSizeFactors(dset)
  dset <- DESeq2::estimateDispersions(dset)
  gene_expr <- DESeq2::getVarianceStabilizedData(dset)
  rownames(gene_expr) <- rownames(raw_counts_integer)

  # Pathways activity
  pathway_activity <- progeny::progeny(gene_expr, scale = FALSE, organism = "Human", verbose = verbose)

  # check what is the percentage of genes we have in our data
  model_pathways <- progeny::getModel(organism = "Human", top = 100)
  full_pathway_sig <- unique(unlist(lapply(colnames(model_pathways), function(pathway) {
    top_genes_pathway <- rownames(model_pathways)[apply(model_pathways, 2, function(X) X != 0)[, pathway]]
    return(top_genes_pathway)
  })))
  genes_kept <- intersect(rownames(gene_expr), full_pathway_sig)
  genes_left <- setdiff(full_pathway_sig, rownames(gene_expr))
  total_genes <- length(genes_left) + length(genes_kept)
  if (verbose) message("Pathway signature genes found in data set: ", length(genes_kept), "/", total_genes, " (", round(length(genes_kept) / total_genes, 3) * 100, "%)")

  if (verbose) message("\nPathway scores computed \n")
  return(as.data.frame(pathway_activity))
}
