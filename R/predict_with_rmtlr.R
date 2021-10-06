#' Predict single-view immune response
#'
#' Obtains single-view predictions of immune response by using a
#' cancer-specific model learned with Regularized
#' Multi-Task Linear Regression algorithm (RMTLR).
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param view_name character string containing the name of the
#' input view.
#' @param view_info character string informing about the family of
#' the input data.
#' @param view_data list containing the data for each input view.
#' @param opt_model_cancer_view_spec cancer-view-specific model
#' feature parameters learned during training.
#' @param opt_xtrain_stats_cancer_view_spec cancer-view-specific
#' features mean and standard deviation of the training set.
#' @param verbose logical flag indicating whether to display
#' messages about the process.
#'
#' @return A lList of predictions matrices, one for each tasks
#' (rows = samples; columns = [runs).
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
#' set.seed(1234)
#' pat_subset <- sample(colnames(RNA_tpm), size = 5)
#' RNA_tpm <- RNA_tpm[, pat_subset]
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' view_name <- "immunecells"
#' view_info <- c(immunecells = "gaussian")
#' view_data <- list(immunecells = as.data.frame(cell_fractions))
#'
#' # Retrieve internal data
#' opt_models <- suppressMessages(easierData::get_opt_models())
#' opt_xtrain_stats <- suppressMessages(easierData::get_opt_xtrain_stats())
#'
#' opt_model_cancer_view_spec <- lapply(view_name, function(X) {
#'     return(opt_models[[cancer_type]][[X]])
#' })
#' names(opt_model_cancer_view_spec) <- view_name
#' opt_xtrain_stats_cancer_view_spec <- lapply(view_name, function(X) {
#'     return(opt_xtrain_stats[[cancer_type]][[X]])
#' })
#' names(opt_xtrain_stats_cancer_view_spec) <- view_name
#'
#' # Predict using rmtlr
#' prediction_view <- predict_with_rmtlr(
#'     view_name = view_name,
#'     view_info = view_info,
#'     view_data = view_data,
#'     opt_model_cancer_view_spec = opt_model_cancer_view_spec,
#'     opt_xtrain_stats_cancer_view_spec = opt_xtrain_stats_cancer_view_spec
#' )
predict_with_rmtlr <- function(view_name,
                               view_info,
                               view_data,
                               opt_model_cancer_view_spec,
                               opt_xtrain_stats_cancer_view_spec,
                               verbose = TRUE) {

    # initialize variables
    P <- length(view_info)
    K <- ncol(opt_model_cancer_view_spec[[1]][[1]])
    standardize_any <- TRUE
    tasks <- names(opt_model_cancer_view_spec[[1]])

    # Algorithm do not deal with NA values: here we removed features with NA values,
    # patients with all NA values are not removed
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
        for (m in seq_len(P)) {

            # Check features availability (we account for mismatches in
            # new releases of progeny, dorothea or quantiseqr)
            keep_pos <- stats::na.omit(match(
                colnames(prediction_X[[m]]),
                rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)
            ))
            keep_names <- intersect(
                colnames(prediction_X[[m]]),
                rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)
            )
            prediction_X[[m]] <- prediction_X[[m]][, keep_names]

            # Normalization should be done taking into account the train set
            opt_xtrain_stats_cancer_view_spec[[m]]$mean <- opt_xtrain_stats_cancer_view_spec[[m]]$mean[keep_pos, ]
            opt_xtrain_stats_cancer_view_spec[[m]]$sd <- opt_xtrain_stats_cancer_view_spec[[m]]$sd[keep_pos, ]

            prediction_X_norm <- lapply(seq_len(K), function(k) {
                prediction_X[[m]] <- calc_z_score(
                    X = prediction_X[[m]], mean = opt_xtrain_stats_cancer_view_spec[[m]]$mean[, k],
                    sd = opt_xtrain_stats_cancer_view_spec[[m]]$sd[, k]
                )
            })
        }
    }

    # perform prediction
    prediction_cv <- lapply(seq_len(K), function(k) {
        coef_matrix <- vapply(tasks, function(task) {
            state[[view_name]][[task]][, k]
        }, FUN.VALUE = numeric(nrow(state[[1]][[1]])))
        rmtlr_test(prediction_X_norm[[k]], coef_matrix)
    })

    predictions_all_tasks_cv <- lapply(tasks, function(task) {
        prediction_task_cv <- vapply(seq_len(K), function(k) {
            prediction_cv[[k]][, task]
        }, FUN.VALUE = numeric(nrow(prediction_X[[1]])))
        return(prediction_task_cv)
    })
    names(predictions_all_tasks_cv) <- tasks
    return(predictions_all_tasks_cv)
}
