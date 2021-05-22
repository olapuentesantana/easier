#' Pre-process data from mariathasan cohort
#'
#' `preprocess_mariathasan` prepares data from mariathasan cohort to be used along the vignette.
#'
#' @param cds CountDataSet object
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return A list with three data.frames:
#' \describe{
#'  \item{counts}{A data.frame containing gene count data with genes as rows and samples as columns.}
#'  \item{tpm}{A data.frame containing tpm data with genes as rows and samples as columns.}
#'  \item{clinical}{A data.frame that includes information on certain clinical variables.}
#' }
#' @export
#'
#' @examples
#'
#' data(cds)
#'
#' data_mariathsan <- preprocess_mariathasan(cds)
#'
preprocess_mariathasan <- function(cds,
                                   verbose = TRUE) {
  # Load CountDataSet object
  data(cds)
  # Retrieve counts and genes information
  gene_count <- counts(cds)
  gene_info <- fData(cds)
  clin_info <- pData(cds)

  # Are there any duplicated symbols?
  all_duplicates <- gene_info[duplicated(gene_info$symbol), ]
  in_duplicates <- all_duplicates[!all_duplicates$symbol %in% c(NA, ""), ]
  out_duplicates <- all_duplicates[all_duplicates$symbol %in% c(NA, ""), ]

  # NCBI source (CSNK1E =1454; PTEP2-CSNK1E = 102800317)
  gene_info[which(gene_info$entrez_id == 102800317), c("symbol", "Symbol")] <- c("TPTEP2-CSNK1E", "TPTEP2-CSNK1E") # NCBI source (CSNK1E =1454; PTEP2-CSNK1E = 102800317)

  # Remove entrez id duplicates (empty symbols)
  gene_info <- gene_info[-which(gene_info$entrez_id %in% out_duplicates$entrez_id), ]
  gene_count <- gene_count[-which(rownames(gene_count) %in% out_duplicates$entrez_id), ]

  # Remove NA symbol values
  where_NA <- is.na(gene_info$symbol)
  gene_info <- gene_info[!where_NA,]
  gene_count <- gene_count[!where_NA,]
  rownames(gene_count) <- gene_info$symbol # Perfect match

  # obtain tpm from counts
  gene_length <- gene_info$length
  gene_count_by_len <- sweep(gene_count, 1, gene_length/1000, FUN = "/")
  scaling_factor <- colSums(gene_count_by_len) / 1e6
  gene_tpm <- sweep(gene_count_by_len, 2, scaling_factor, FUN = "/")

  # keep only those patients with available and unambiguous response (CR and PD)
  BOR <- clin_info[, "Best Confirmed Overall Response", drop = FALSE]
  BOR_subset <- BOR[!BOR$`Best Confirmed Overall Response` == "NE", , drop = FALSE]
  BOR_subset <- BOR_subset[BOR_subset$`Best Confirmed Overall Response` %in% c("CR", "PD"), , drop = FALSE]
  gene_tpm_subset <- as.data.frame(gene_tpm[, match(rownames(BOR_subset), colnames(gene_tpm))])
  gene_count_subset <- as.data.frame(gene_count[, match(rownames(BOR_subset), colnames(gene_count))])
  clin_info_subset <- as.data.frame(clin_info[rownames(clin_info) %in% rownames(BOR_subset), , drop = FALSE])

  data <- list(counts=gene_count_subset, tpm=gene_tpm_subset, clinical=clin_info_subset)

  if (verbose) message("Mariathasan data retrieved successfully!")
  return(data)
}
