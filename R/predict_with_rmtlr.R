#' Collects Regularized Multi-Task Linear Regression (RMTLR) cancer-specfic model predictions
#'
#' This function puts together predicted immune response by using a cancer-specifc model learned
#' from training with RMTLR algorithm.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param view_name character string containing the name of the input view.
#' @param view_info character string informing about the family of the input data.
#' @param view_data list containing the data for each input view.
#' @param cancer_type character string indicating cancer type to specify cancer-specific optimization model
#' to be used.
#' @param verbose logical flag indicating whether to display messages about the process.
#'
#' @return A lList of predictions matrices, one for each tasks (rows = samples; columns = runs).
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = gene_tpm)
#'
#' view_name <- "immunecells"
#' view_info <- c(immunecells = "gaussian")
#' view_data <- list(immunecells = as.data.frame(cell_fractions))
#'
#' # Predict using rmtlr
#' prediction_view <- predict_with_rmtlr(
#'   view_name = view_name,
#'   view_info = view_info,
#'   view_data = view_data,
#'   cancer_type = "SKCM"
#' )
predict_with_rmtlr <- function(view_name,
                               view_info,
                               view_data,
                               cancer_type,
                               verbose = TRUE) {

  # Initialize variables
  opt_model_cancer_view_spec <- lapply(view_name, function(X) {
    return(opt_models[[cancer_type]][[X]])
  })
  names(opt_model_cancer_view_spec) <- view_name
  opt_xtrain_stats_cancer_view_spec <- lapply(view_name, function(X) {
    return(opt_xtrain_stats[[cancer_type]][[X]])
  })
  names(opt_xtrain_stats_cancer_view_spec) <- view_name

  P <- length(view_info)
  K <- ncol(opt_model_cancer_view_spec[[1]][[1]])
  standardize_any <- TRUE
  tasks <- names(opt_model_cancer_view_spec[[1]])

  # Algorithm do not deal with NA values: here we removed features with NA values, patients with all NA values are not removed
  view_data_new <- lapply(names(view_data), function(x) {
    tmp_data <- view_data[[x]]
    if (all(is.na(apply(tmp_data, 1, sum)))) {
      NA_sum <- apply(tmp_data, 2, sum)
      tmp_data <- tmp_data[, !is.na(NA_sum)]
    }
    return(tmp_data)
  })
  names(view_data_new) <- names(view_data)

  state <- opt_model_cancer_view_spec
  prediction_X <- view_data_new

  # standardize
  if (standardize_any == TRUE) {
    for (m in 1:P) {

      # Check features availability
      keep_pos <- stats::na.omit(match(colnames(prediction_X[[m]]), rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)))
      keep_names <- intersect(colnames(prediction_X[[m]]), rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean))
      prediction_X[[m]] <- prediction_X[[m]][, keep_names]

      # Normalization should be done taking into account the train set
      opt_xtrain_stats_cancer_view_spec[[m]]$mean <- opt_xtrain_stats_cancer_view_spec[[m]]$mean[keep_pos, ]
      opt_xtrain_stats_cancer_view_spec[[m]]$sd <- opt_xtrain_stats_cancer_view_spec[[m]]$sd[keep_pos, ]

      prediction_X_norm <- lapply(1:K, function(k) {
        prediction_X[[m]] <- calc_z_score(
          X = prediction_X[[m]], mean = opt_xtrain_stats_cancer_view_spec[[m]]$mean[, k],
          sd = opt_xtrain_stats_cancer_view_spec[[m]]$sd[, k]
        )
      })
    }
  }

  # perform prediction
  prediction_cv <- lapply(1:K, function(k) {
    coef_matrix <- sapply(tasks, function(task) {
      state[[view_name]][[task]][, k]
    })
    rmtlr_test(prediction_X_norm[[k]], coef_matrix)
  })

  # save predictions
  prediction_cv <- lapply(1:K, function(k) {
    coef_matrix <- sapply(tasks, function(task) {
      state[[view_name]][[task]][, k]
    })
    rmtlr_test(prediction_X_norm[[k]], coef_matrix)
  })

  predictions_all_tasks_cv <- lapply(tasks, function(task) {
    prediction_task_cv <- sapply(1:K, function(k) {
      prediction_cv[[k]][, task]
    })
    return(prediction_task_cv)
  })
  names(predictions_all_tasks_cv) <- tasks
  return(predictions_all_tasks_cv)
}
