#' Immune response prediction
#'
#' `predict_immune_response` predicts immune response using two algorithms:
#' regularized multi-task linear regression (RMTLR) and bayesian efficient multi-kernel algorithm (BEMKL). While
#' BEMKL can exploit information across different input and output datasets,
#' RMTLR can only do so for response variables.
#' Another advantage of BEMKL is missing data handling, which is not the case
#' for the other algorithm.
#'
#' These algorithms use model parameters learned during training on different
#' types of data in order to compute the immune response.
#'
#' @importFrom utils combn
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param pathways numeric matrix with data
#' @param immunecells numeric matrix with data
#' @param tfs numeric matrix with data
#' @param lrpairs numeric matrix with data
#' @param ccpairs numeric matrix with data
#' @param cancertype string character
#'
#' @return Predictions for each model building.
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = Riaz_data$tpm_RNAseq)
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathways_scores(RNA_counts = Riaz_data$raw_counts_RNAseq,
#' remove_genes_ICB_proxies = TRUE)
#'
#' # Computation of TF activity
#' tf_activity <- compute_TF_activity(RNA_tpm = Riaz_data$tpm_RNAseq,
#' remove_genes_ICB_proxies = FALSE)
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(RNA_tpm = Riaz_data$tpm_RNAseq,
#' remove_genes_ICB_proxies = FALSE, cancertype = "pancan")
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs_grouped(lrpairs = lrpairs_weights$LRpairs,
#' cancertype = "pancan")
#'
#' # Predict patients' immune response
#' predictions_immune_response <- predict_immune_response(pathways = pathways_activity,
#' immunecells = cell_fractions,
#' lrpairs = lrpairs_weights,
#' tfs = tf_activity
#' ccpairs = ccpairs_scores,
#' cancertype = "SKCM")
predict_immune_response <- function(pathways = NULL,
                                    immunecells = NULL,
                                    tfs = NULL,
                                    lrpairs = NULL,
                                    ccpairs = NULL,
                                    cancertype) {
  if (missing(cancertype)) stop("cancer type needs to be specified")
  if (all(is.null(pathways), is.null(immunecells), is.null(tfs), is.null(lrpairs), is.null(ccpairs))) stop("none signature specified")

  # Simplify efforts: get data in lowercase variables
  pathways.cor <- pathways
  lrpairs.spec.pc <- lrpairs
  ccpairsgroupedscores.spec.pc <- ccpairs

  # Initialize variables
  views <- c(
    Pathways.cor = "gaussian",
    ImmuneCells = "gaussian",
    TFs = "gaussian",
    LRpairs.spec.pc = "gaussian",
    CCpairsGroupedScores.spec.pc = "gaussian"
  )

  view_combinations <- NULL

  algorithm <- c("RMTLR") # ,"BEMKL")

  # Check which views are missing
  miss_views <- c(
    ifelse(missing(pathways), NA, 1),
    ifelse(missing(immunecells), NA, 2),
    ifelse(missing(tfs), NA, 3),
    ifelse(missing(lrpairs), NA, 4),
    ifelse(missing(ccpairs), NA, 5)
  )

  # Possible combinations
  possible_combo <- combn(miss_views, m = 2)[, 1:9]

  # Remove combinations with are not feasible due to missing views
  if (anyNA(miss_views)) {
    possible_combo <- possible_combo[, !is.na(colSums(possible_combo)), drop = FALSE]
  }

  # Views single
  view_simples <- lapply(miss_views[!is.na(miss_views)], function(X) {
    tmp <- views[X]
    return(tmp)
  })

  # Views combination
  if (is.matrix(possible_combo) & dim(possible_combo)[2] > 1) {
    view_combinations <- lapply(1:ncol(possible_combo), function(X) {
      tmp <- views[possible_combo[, X]]
      return(tmp)
    })
  }
  view_combinations <- c(view_simples, view_combinations)

  # Remove unavailable combo
  combo_names <- sapply(1:length(view_combinations), function(X) {
    paste(names(view_combinations[[X]]), collapse = "_")
  })

  # Immune cells features curation:
  if (missing(immunecells) == FALSE) {
    colnames(immunecells) <- gsub(".", "_", colnames(immunecells), fixed = TRUE)
  }

  all_predictions <- lapply(1:length(view_combinations), function(X) {
    view_info <- view_combinations[[X]]
    view_name <- paste(names(view_info), collapse = "_")
    view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
    names(view_data) <- names(view_info)
    message(X, ".view source: ", view_name, "\n")

    # Predict immune response using model parameters
    summary_alg <- lapply(algorithm, function(alg) {
      if (alg %in% c("BEMKL")) {
        pred_alg <- predict_with_bemkl(
          view_name = view_name,
          view_info = view_info,
          view_data = view_data,
          learned_model = trained_models[[cancertype]][[view_name]]
        )
      } else if (alg %in% c("RMTLR")) {
        pred_alg <- predict_with_rmtlr(
          view_name = view_name,
          view_info = view_info,
          view_data = view_data,
          learned_model = trained_models[[cancertype]][[view_name]]
        )
      }
      return(pred_alg)
    })
    names(summary_alg) <- algorithm
    return(summary_alg)
  })
  names(all_predictions) <- combo_names
  return(all_predictions)
}
