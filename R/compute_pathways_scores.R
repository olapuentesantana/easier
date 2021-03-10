#' Compute pathway activity from gene expression using PROGENy
#'
#' This function infers pathway activity from gene expression in raw counts from bulk RNA-seq data
#' using PROGENy method from (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018).
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#' getVarianceStabilizedData rlog
#' @import progeny
#' @importFrom stats na.exclude
#'
#' @param RNA_counts A data.frame containing raw counts values with HGNC symbols in rows and samples in columns.
#' @param remove_genes_ICB_proxies A logical value indicating whether to remove signature genes involved
#' in the derivation of hallmarks of immune response.
#' @param verbose A logical value indicating whether to display messages about the number of pathway signature
#' genes found in the gene expression data provided.
#'
#' @return A matrix with samples in rows and pathways in columns.
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#' library("progeny")
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathways_scores(
#'   RNA_counts = Riaz_data$raw_counts_RNAseq,
#'   remove_genes_ICB_proxies = TRUE)
#' head(pathway_activity)
compute_pathways_scores <- function(RNA_counts,
                                    remove_genes_ICB_proxies = TRUE,
                                    verbose = TRUE) {

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
    raw_counts.integer <- apply(raw_counts, 2, as.integer)
    rownames(raw_counts.integer) <- rownames(raw_counts)
  } else {
    raw_counts.integer <- raw_counts
  }

  # Variance stabilizing transformation (DESeq2 package)
  # Integer count matrix, a data frame with the sample info,  design =~1 to consider all samples as part of the same group.

  # Column data:
  colData <- data.frame(id = colnames(raw_counts.integer))

  if (verbose) message("DESeq2 normalization -->\n")
  # Construction a DESeqDataSet: (Forced all to be data.frames($ operator))
  dset <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts.integer,
    colData = colData,
    design = ~1
  )

  # Variance stabilization transformation
  dset <- DESeq2::estimateSizeFactors(dset)
  dset <- DESeq2::estimateDispersions(dset)
  # gene_expr <- DESeq2::rlog(raw_counts.integer)
  gene_expr <- DESeq2::getVarianceStabilizedData(dset)
  rownames(gene_expr) <- rownames(raw_counts.integer)

  # Pathways activity (Progeny package)
  # library(progeny)
  pathway_activity <- progeny::progeny(gene_expr, scale = FALSE, organism = "Human", verbose = verbose)

  # check what is the percentage of genes we have in our data
  all_pathway_responsive_genes <- unique(unlist(top_100_per_pathway_responsive_genes))
  genes_kept <- intersect(rownames(gene_expr), all_pathway_responsive_genes)
  genes_left <- setdiff(all_pathway_responsive_genes, rownames(gene_expr))
  total_genes <- length(genes_left) + length(genes_kept)
  if (verbose) message("Pathway signature genes found in data set: ", length(genes_kept), "/", total_genes, " (", round(length(genes_kept)/total_genes, 3) * 100,"%)")

  if (verbose) message("\n Pathway scores computed \n")
  return(as.data.frame(pathway_activity))

}
