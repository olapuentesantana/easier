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
#' @importFrom stats na.exclude
#' @importFrom dplyr filter
#' @importFrom progeny progeny getModel
#' @importFrom decoupleR get_progeny run_wmean
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom easierData get_cor_scores_genes
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC
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
compute_pathway_activity <- function(RNA_tpm = NULL,
                                     remove_sig_genes_immune_response = TRUE,
                                     verbose = TRUE) {
  # Some checks
  if (is.null(RNA_tpm)) stop("TPM gene expression data not found")
  
  # Retrieve internal data
  cor_scores_genes <- suppressMessages(easierData::get_cor_scores_genes())

  # Gene expression data
  tpm <- RNA_tpm
  genes <- rownames(tpm)
  
  # HGNC symbols are required
  if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE)
  
  gene_expr <- t(tpm)
  # redefine gene names to match TF-target network
  E <- t(gene_expr)
  newNames <- gsub(".", "-", rownames(E), fixed = TRUE)
  rownames(E) <- newNames
  
  ## progeny network
  net <- decoupleR::get_progeny(organism = 'human', top = 100)

  all_regulated_transcripts <- unique(net$target)

  # check what is the percentage of genes we have in our data
  genes_kept <- intersect(rownames(E), all_regulated_transcripts)
  genes_left <- setdiff(all_regulated_transcripts, rownames(E))
  
  # check what is the percentage of regulated transcripts that we have in our data
  if (verbose) {
    message(
      "Regulated transcripts found in data set: ", length(genes_kept), "/",
      length(all_regulated_transcripts), " (",
      round(length(genes_kept) / length(all_regulated_transcripts), 3) * 100, "%)"
    )
  }
  
  # Remove list of genes used to build proxy's of ICB response
  if (remove_sig_genes_immune_response) {
    if (verbose) message("Removing signature genes of hallmarks of immune response \n")
    idy <- stats::na.exclude(match(cor_scores_genes, rownames(E)))
    E <- E[-idy, ]
  }

  # Remove genes with all NA/Inf values
  E <- E[!is.na(apply(E, 1, sum)), ]
  E <- E[!is.infinite(apply(E, 1, sum)), ]
  
  # Pathway activity: run wmean
  pathway_activity_df <- decoupleR::run_wmean(mat = E,
                                              net = net,
                                              .source='source',
                                              .target='target',
                                              .mor='weight', 
                                              minsize = 5)
  # To matrix
  pathway_activity <- pathway_activity_df %>%
    dplyr::filter(statistic == "wmean") %>%
    tidyr::pivot_wider(id_cols = condition,
                       names_from = source,
                       values_from = score) %>%
    tibble::column_to_rownames("condition")

  if (verbose) message("\nPathway activity scores computed! \n")
  return(as.data.frame(pathway_activity))
}
