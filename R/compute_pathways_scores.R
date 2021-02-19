#' Compute pathways score
#'
#' `compute_pathways_scores` infers pathway activity from raw counts RNA-seq data.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#' getVarianceStabilizedData rlog
#' @import progeny
#' @importFrom stats na.exclude
#'
#' @param RNA_counts numeric matrix of read counts with rows=genes and columns=samples
#' @param remove_genes_ICB_proxies boolean variable to remove all those genes involved
#' in the computation of ICB proxy's of response
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{scores}{pathway activity matrix with rows=samples and columns=pathways}
#'     \item{transcripts_kept}{vector with available gene names}
#'     \item{transcripts_left}{vector with missing gene names}
#'   }
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathways_scores(RNA.counts = Riaz_data$raw_counts_RNAseq,
#' remove.genes.ICB_proxies = TRUE)
compute_pathways_scores <- function(RNA_counts,
                                    remove_genes_ICB_proxies = TRUE) {

  # Gene expression data
  raw_counts <- RNA_counts
  genes <- rownames(raw_counts)

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE))

  # Remove list of genes used to build proxy's of ICB response
  if (remove_genes_ICB_proxies) {
    message("Removing signatures genes for proxy's of ICB response  \n")
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

  message("DESeq2 -->  \n")
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
  Pathway_scores <- progeny::progeny(gene_expr, scale = FALSE, organism = "Human", verbose = TRUE)

  # check what is the percentage of genes we have in our data
  all_pathway_responsive_genes <- unique(unlist(top_100_per_pathway_responsive_genes))
  genes_kept <- intersect(rownames(gene_expr), all_pathway_responsive_genes)
  genes_left <- setdiff(all_pathway_responsive_genes, rownames(gene_expr))

  # Output list:
  Pathways <- list(
    scores = as.data.frame(Pathway_scores),
    transcripts_kept = length(genes_kept),
    transcripts_left = length(genes_left)
  )

  message("\n Pathway scores computed \n")
  return(Pathways)
}
