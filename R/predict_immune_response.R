#' Function to compute predicted immune response
#'
#' `predict_immune_response` predicts patients' immune response using regularized multi-task linear regression (RMTLR).
#' The predicted immune response is calculated based on the input features data and the parameters learned during model training.
#'
#' @importFrom utils combn
#' @importFrom stats na.omit
#' @importFrom BiocParallel register bplapply MulticoreParam
#'
#' @export
#'
#' @param pathways A numeric matrix (rows = samples; columns = pathways).
#' @param immunecells A numeric matrix (rows = samples; columns = cell types).
#' @param tfs A numeric matrix (rows = samples; columns = transcription factors).
#' @param lrpairs A numeric matrix (rows = samples; columns = ligand-receptor pairs).
#' @param ccpairs A numeric matrix (rows = samples; columns = cell-cell pairs).
#' @param cancer_type A character string indicating which cancer-specific model should be used to compute the predictions.
#' @param verbose A logical flag indicating whether to display messages about the process.
#'
#' @return A list containing the predictions for each quantitative descriptor and for each task.
#' Given that the model training was repeated 100 times with randomized-cross validation, a set of 100 predictions is returned.
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_count <- mariathasan_data$counts
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Computation of cell fractions (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = gene_tpm)
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019; Schubert et al., Nat Commun, 2018)
#' pathway_activity <- compute_pathways_scores(RNA_counts = gene_count,
#' remove_genes_ICB_proxies = TRUE)
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activity <- compute_TF_activity(RNA_tpm = gene_tpm,
#' remove_genes_ICB_proxies = FALSE)
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(RNA_tpm = gene_tpm,
#' remove_genes_ICB_proxies = FALSE,
#' cancer_type = "pancan")
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs_grouped(lrpairs = lrpair_weights,
#' cancer_type = "pancan")
#'
#' # Predict patients' immune response
#' predictions_immune_response <- predict_immune_response(pathways = pathway_activity,
#' immunecells = cell_fractions,
#' tfs = tf_activity,
#' lrpairs = lrpair_weights,
#' ccpairs = ccpair_scores,
#' cancer_type = "BLCA")
predict_immune_response <- function(pathways = NULL,
                                    immunecells = NULL,
                                    tfs = NULL,
                                    lrpairs = NULL,
                                    ccpairs = NULL,
                                    cancer_type,
                                    verbose = TRUE) {
  if (missing(cancer_type)) stop("cancer type needs to be specified")
  if (all(is.null(pathways), is.null(immunecells), is.null(tfs), is.null(lrpairs), is.null(ccpairs))) stop("none signature specified")

  # Initialize variables
  views <- c(
    pathways = "gaussian",
    immunecells = "gaussian",
    tfs = "gaussian",
    lrpairs = "gaussian",
    ccpairs = "gaussian"
  )

  # Check which views are missing
  miss_views <- c(
    ifelse(missing(pathways), NA, 1),
    ifelse(missing(immunecells), NA, 2),
    ifelse(missing(tfs), NA, 3),
    ifelse(missing(lrpairs), NA, 4),
    ifelse(missing(ccpairs), NA, 5)
  )

  # Single views
  view_simples <- lapply(miss_views[!is.na(miss_views)], function(X) {
    tmp <- views[X]
    return(tmp)
  })

  # All corresponding views
  view_combinations <- view_simples

  compute_prediction <- function(view, verbose){

    view_info <- view_combinations[[view]]
    view_name <- paste(names(view_info), collapse = "_")
    view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
    names(view_data) <- names(view_info)
    if (verbose) message("Computing predictions using ", view_name, "...\n")

    # Predict immune response using RMTLR model parameters
    prediction_view <- predict_with_rmtlr(
      view_name = view_name,
      view_info = view_info,
      view_data = view_data,
      cancer_type = cancer_type
    )

    return(prediction_view)
  }
  # Parallelize views model predictions
  BiocParallel::register(BiocParallel::MulticoreParam(workers = 2))
  all_predictions <- BiocParallel::bplapply(1:length(view_combinations), FUN = compute_prediction, verbose = verbose)

  names(all_predictions) <- sapply(1:length(view_combinations), function(X) {
    paste(names(view_combinations[[X]]), collapse = "_")
  })
  return(all_predictions)
}
