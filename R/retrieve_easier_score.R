#' Retrieve easier scores of immune response
#'
#' Calculates easier score and if applicable, both weighted average and
#' penalized score based on the combination of easier score and TMB.
#'
#' @param predictions_immune_response list containing the predictions
#' for each quantitative descriptor and for each task. This is the
#' output from \code{predict_immune_response}.
#' @param TMB_values numeric vector containing patients' tumor
#' mutational burden (TMB) values.
#' @param easier_with_TMB character string indicating which approach
#' should be used to integrate easier with TMB: "weighted_average"
#' (default) and "penalized_score".
#' @param weight_penalty integer value from 0 to 1, which is used to
#' define the weight or penalty for combining easier and TMB scores
#' based on a weighted average or penalized score, in order to derive
#' a score of patient's likelihood of immune response. The default
#' value is 0.5.
#' @param verbose logical flag indicating whether to display messages
#' about the process.
#'
#' @return A data.frame with samples in rows and easier scores in columns.
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
#' cancer_type <- metadata(dataset_mariathasan)[["cancertype"]]
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'   "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'   "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activities <- compute_TF_activity(
#'   RNA_tpm = RNA_tpm
#' )
#'
#' # Predict patients' immune response
#' predictions <- predict_immune_response(
#'   tfs = tf_activities,
#'   cancer_type = cancer_type,
#'   verbose = TRUE
#' )
#'
#' # retrieve clinical response
#' patient_ICBresponse <- colData(dataset_mariathasan)[["BOR"]]
#' names(patient_ICBresponse) <- colData(dataset_mariathasan)[["pat_id"]]
#'
#' # retrieve TMB
#' TMB <- colData(dataset_mariathasan)[["TMB"]]
#' names(TMB) <- colData(dataset_mariathasan)[["pat_id"]]
#'
#' patient_ICBresponse <- patient_ICBresponse[names(patient_ICBresponse) %in% pat_subset]
#' TMB <- TMB[names(TMB) %in% pat_subset]
#'
#' easier_derived_scores <- retrieve_easier_score(
#'   predictions_immune_response = predictions,
#'   TMB_values = TMB,
#'   easier_with_TMB = c("weighted_average", "penalized_score"),
#'   weight_penalty = 0.5
#' )
#' \donttest{
#'
#' RNA_counts <- assays(dataset_mariathasan)[["counts"]]
#' RNA_counts <- RNA_counts[, colnames(RNA_counts) %in% pat_subset]
#'
#' # Computation of cell fractions  (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019;
#' # Schubert et al., Nat Commun, 2018)
#' pathway_activities <- compute_pathway_activity(
#'   RNA_counts = RNA_counts,
#'   remove_sig_genes_immune_response = TRUE
#' )
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(
#'   RNA_tpm = RNA_tpm,
#'   cancer_type = "pancan"
#' )
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs(
#'   lrpairs = lrpair_weights,
#'   cancer_type = "pancan"
#' )
#'
#' # Predict patients' immune response
#' predictions <- predict_immune_response(
#'   pathways = pathway_activities,
#'   immunecells = cell_fractions,
#'   tfs = tf_activities,
#'   lrpairs = lrpair_weights,
#'   ccpairs = ccpair_scores,
#'   cancer_type = cancer_type,
#'   verbose = TRUE
#' )
#'
#' easier_derived_scores <- retrieve_easier_score(
#'   predictions_immune_response = predictions,
#'   TMB_values = TMB,
#'   easier_with_TMB = c("weighted_average", "penalized_score"),
#'   weight_penalty = 0.5
#' )
#' }
retrieve_easier_score <- function(predictions_immune_response = NULL,
                                  TMB_values = NULL,
                                  easier_with_TMB = c("weighted_average", "penalized_score"),
                                  weight_penalty,
                                  verbose = TRUE) {

  # Some checks
  if (is.null(predictions_immune_response)) stop("None predictions found")
  if (is.null(TMB_values)) {
    TMB_available <- FALSE
    easier_with_TMB <- "none"
    patients_to_keep <- rownames(predictions_immune_response[[1]][[1]])
  } else {
    TMB_available <- TRUE
    if (missing(easier_with_TMB)) {
      easier_with_TMB <- "weighted_average"
    }
    if (anyNA(TMB_values)) {
      warning(
        "TMB data contains NA values, ",
        "excluding those patients.."
      )
    }
    patients_to_keep <- names(TMB_values[!is.na(TMB_values)])
    TMB_values <- as.numeric(TMB_values[patients_to_keep])
    names(TMB_values) <- patients_to_keep
  }

  # Initialize function
  if ("weighted_average" %in% easier_with_TMB) {
    linear_func <- function(x) {
      min_x <- min(x)
      max_x <- max(x)
      x01 <- (x - min_x) / (max_x - min_x)
      return(x01)
    }
  }

  # Initialize variables
  # view names
  views <- names(predictions_immune_response)
  # task names
  tasks <- names(predictions_immune_response[[1]])

  # Ensemble predictor
  ## across runs
  ensemble_df <- lapply(views, function(spec_view) {
    ensemble_df <- vapply(tasks, function(spec_task) {
      df <- predictions_immune_response[[spec_view]][[spec_task]]
      df <- df[patients_to_keep, ]
      df_runs <- rowMeans(df)
    }, FUN.VALUE = numeric(length(patients_to_keep)))
    return(ensemble_df)
  })
  names(ensemble_df) <- views
  ## across views
  overall_df <- vapply(tasks, function(spec_task) {
    overall_df <- apply(cbind(
      ensemble_df$pathways[, spec_task],
      ensemble_df$immunecells[, spec_task],
      ensemble_df$tfs[, spec_task],
      ensemble_df$lrpairs[, spec_task],
      ensemble_df$ccpairs[, spec_task]
    ), 1, mean)
  }, FUN.VALUE = numeric(length(patients_to_keep)))
  easier_score <- apply(overall_df, 1, mean)

  # save into data.frame
  scores <- data.frame(easier_score = easier_score)

  if (!"none" %in% easier_with_TMB) {
    rp_df <- data.frame(
      prediction_easier = easier_score,
      patient = rownames(overall_df)
    )
    ### default weight_penalty value = 0.5
    if (missing(weight_penalty)) weight_penalty <- 0.5
    ### add TMB
    rp_df$TMB <- TMB_values
    ### categorize TMB
    if (length(unique(rp_df$TMB)) > 3) {
      rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
      rp_df$TMB <- rp_df$TMBcat
    }
    ### compute the integrated score as weighted average
    if ("weighted_average" %in% easier_with_TMB) {
      pred_lin <- linear_func(rp_df$prediction_easier)
      names(pred_lin) <- rp_df$patient
      TMB_lin <- linear_func(rp_df$TMB)
      names(TMB_lin) <- rp_df$patient
      pred_averaged_rf <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_averaged <- apply(cbind((1 - p) * pred_lin, (p) * TMB_lin), 1, mean)
      }, FUN.VALUE = numeric(length(pred_lin)))
      weight_penalty_pos <- match(weight_penalty, (seq_len(11) - 1) / 10)
      w_avg_score <- pred_averaged_rf[, weight_penalty_pos]
      # save into data.frame
      scores <- cbind(scores, data.frame(w_avg_score = w_avg_score))
    }
    ### compute the integrated score for different penalties
    if ("penalized_score" %in% easier_with_TMB) {
      pred_combined <- rp_df$prediction_easier
      names(pred_combined) <- rp_df$patient
      pred_combined_rf <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
        pred_combined[rp_df$TMBcat == 1] <- pred_combined[rp_df$TMBcat == 1] - p
        pred_combined[rp_df$TMBcat == 3] <- pred_combined[rp_df$TMBcat == 3] + p
        return(pred_combined)
      }, FUN.VALUE = numeric(length(pred_combined)))
      weight_penalty_pos <- match(weight_penalty, (seq_len(11) - 1) / 10)
      pen_score <- pred_combined_rf[, weight_penalty_pos]
      # save into data.frame
      scores <- cbind(scores, data.frame(pen_score = pen_score))
    }
  }
  return(scores)
}
