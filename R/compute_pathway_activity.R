#' Compute pathway activity from gene expression using PROGENy
#'
#' Infers pathway activity from counts bulk gene expression
#' using PROGENy method from Holland et al., BBAGRM, 2019 and
#' Schubert et al., Nat Commun, 2018.
#'
#' @references Schubert M, Klinger B, Klunemann M, Sieber A, Uhlitz F, Sauer S,
#' Garnett MJ, Bluthgen N, Saez-Rodriguez J. “Perturbation-response genes reveal
#' signaling footprints in cancer gene expression.” Nature Communications:
#' 10.1038/s41467-017-02391-6
#'
#' @references Holland CH, Szalai B, Saez-Rodriguez J. "Transfer of regulatory
#' knowledge from human to mouse for functional genomics analysis." Biochimica et
#' Biophysica Acta (BBA) - Gene Regulatory Mechanisms. 2019.
#' DOI: 10.1016/j.bbagrm.2019.194431.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#' getVarianceStabilizedData
#' @importFrom stats na.exclude
#' @importFrom progeny progeny getModel
#' @importFrom easierData get_cor_scores_genes
#'
#' @param RNA_counts data.frame containing raw counts values with HGNC
#' gene symbols as row names and samples identifiers as column names.
#' @param remove_sig_genes_immune_response logical value indicating
#' whether to remove signature genes involved in the derivation of
#' hallmarks of immune response. This list is available from easierData
#' package through \code{easierData::get_cor_scores_genes()}.
#' @param verbose logical value indicating whether to display messages
#' about the number of pathway signature genes found in the gene
#' expression data provided.
#'
#' @return A matrix of activity scores with samples in rows and pathways
#' in columns.
#'
#' @export
#'
#' @examples
#' # using a SummarizedExperiment object
#' library(SummarizedExperiment)
#' # Using example exemplary dataset (Mariathasan et al., Nature, 2018)
#' # from easierData. Original processed data is available from
#' # IMvigor210CoreBiologies package.
#' library("easierData")
#'
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_counts <- assays(dataset_mariathasan)[["counts"]]
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'   "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'   "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_counts <- RNA_counts[, colnames(RNA_counts) %in% pat_subset]
#'
#' # Computation of pathway activity
#' # (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018)
#' pathway_activity <- compute_pathway_activity(
#'   RNA_counts = RNA_counts,
#'   remove_sig_genes_immune_response = TRUE
#' )
compute_pathway_activity <- function(RNA_counts = NULL,
                                     remove_sig_genes_immune_response = TRUE,
                                     verbose = TRUE) {
  # Some checks
  if (is.null(RNA_counts)) stop("Counts gene expression data not found")

  # Retrieve internal data
  cor_scores_genes <- suppressMessages(easierData::get_cor_scores_genes())

  # Gene expression data
  raw_counts <- RNA_counts
  genes <- rownames(raw_counts)

  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) {
    stop("Hgnc gene symbols are required",
      call. = FALSE
    )
  }

  # Remove list of genes used to build proxy's of ICB response
  if (remove_sig_genes_immune_response) {
    if (verbose) message("Removing signature genes of hallmarks of immune response \n")
    idy <- stats::na.exclude(match(cor_scores_genes, rownames(raw_counts)))
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
  # Integer count matrix, a data frame with the sample info,
  # design =~1 to consider all samples as part of the same group.

  # Column data:
  colData <- data.frame(id = colnames(raw_counts_integer))

  if (verbose) message("Gene counts normalization with DESeq2:")
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
  pathway_activity <- progeny::progeny(gene_expr,
    scale = FALSE,
    organism = "Human", verbose = verbose
  )

  # check what is the percentage of genes we have in our data
  model_pathways <- progeny::getModel(organism = "Human", top = 100)
  full_pathway_sig <- unique(unlist(lapply(
    colnames(model_pathways),
    function(pathway) {
      top_genes_pathway <- rownames(model_pathways)[apply(
        model_pathways, 2,
        function(X) X != 0
      )[, pathway]]
      return(top_genes_pathway)
    }
  )))
  genes_kept <- intersect(rownames(gene_expr), full_pathway_sig)
  genes_left <- setdiff(full_pathway_sig, rownames(gene_expr))
  total_genes <- length(genes_left) + length(genes_kept)
  if (verbose) {
    message(
      "Pathway signature genes found in data set: ",
      length(genes_kept), "/",
      total_genes, " (",
      round(length(genes_kept) / total_genes, 3) * 100,
      "%)"
    )
  }

  if (verbose) message("\nPathway activity scores computed! \n")
  return(as.data.frame(pathway_activity))
}
