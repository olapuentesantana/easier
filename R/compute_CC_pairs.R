#' Compute cell-cell interactions scores using computed
#' ligand-receptor weights
#'
#' Infers scores of cell-cell interactions in the tumor
#' microenvironment (Lapuente-Santana et al., Patterns, 2021) using
#' the ligand-receptor weights obtained from \code{compute_LR_pairs}
#' as input.
#'
#' @references Oscar Lapuente-Santana, Maisa van Genderen,
#' Peter A. J. Hilbers, Francesca Finotello, and Federica Eduati.
#' 2021. Interpretable Systems Biomarkers Predict Response to
#' Immune-Checkpoint Inhibitors. Patterns, 100293.
#' https://doi.org/10.1016/j.patter.2021.100293.
#'
#' @importFrom easierData get_lr_frequency_TCGA get_intercell_networks
#'
#' @param lrpairs output of the compute_LR_pairs function. A matrix
#' of log2(TPM +1) weights with samples in rows and ligand-receptor
#' pairs in columns. This is the output from \code{compute_LR_pairs}.
#' @param cancer_type string detailing the cancer type whose cell-cell
#' interaction network will be used. By default, a pan-cancer network
#' is selected whose network represents the union of all ligand-receptor
#' pairs present across the 18 cancer types studied in
#' Lapuente-Santana et al., Patterns, 2021.
#' @param verbose logical value indicating whether to display
#' informative messages about the process.
#'
#' @return A matrix of scores with samples in rows and cell-cell
#' pairs in columns.
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
#' RNA_tpm <- assays(dataset_mariathasan)[["tpm"]]
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'     "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'     "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(
#'     RNA_tpm = RNA_tpm,
#'     cancer_type = "pancan"
#' )
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs(
#'     lrpairs = lrpair_weights,
#'     cancer_type = "pancan"
#' )
compute_CC_pairs <- function(lrpairs = NULL,
                             cancer_type = "pancan",
                             verbose = TRUE) {
    # Some checks
    if (is.null(lrpairs)) stop("ligand-receptor pair weights not found")

    # Retrieve internal data
    lr_frequency <- suppressMessages(easierData::get_lr_frequency_TCGA())
    intercell_networks <- suppressMessages(easierData::get_intercell_networks())

    # remove ligand receptor pairs that are always NA
    na_lrpairs <- apply(lrpairs, 2, function(x) {
        all(is.na(x))
    })
    lrpairs <- lrpairs[, na_lrpairs == FALSE]

    # Binarize the data: set a threshold to 10 TPM
    # Only pairs where both ligand and receptor have TPM > 10 are kept
    lrpairs_binary <- ifelse(lrpairs > log2(10 + 1), 1, 0)

    # keep only the LR.pairs for which I have (non-zero) frequencies in the TCGA
    lrpairs_binary <- lrpairs_binary[, colnames(lrpairs_binary) %in% names(lr_frequency)]

    # cancer type specific network
    intercell_network <- intercell_networks[[cancer_type]]

    # compute the CC score for each patient
    celltypes <- unique(c(
        as.character(intercell_network$cell1),
        as.character(intercell_network$cell2)
    ))

    CC_pairs_score <- do.call(cbind, lapply(celltypes, function(celltype1) {
        do.call(cbind, lapply(celltypes, function(celltype2) {
            compute_CCpair_score(celltype1, celltype2, intercell_network,
                lrpairs_binary, lr_frequency,
                compute_log = TRUE
            )
        }))
    }))

    metadata_CC_pairs <- do.call(rbind, lapply(celltypes, function(celltype1) {
        do.call(rbind, lapply(celltypes, function(celltype2) {
            data.frame(
                CCpair = gsub(" ", "", paste(celltype1, celltype2, sep = "_")),
                celltype1 = celltype1, celltype2 = celltype2
            )
        }))
    }))

    colnames(CC_pairs_score) <- metadata_CC_pairs$CCpair

    if (verbose) message("CC pairs computed \n")
    return(as.data.frame(CC_pairs_score))
}
