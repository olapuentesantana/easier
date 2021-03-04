#' Make predictions using RMTLR
#'
#' `predict_with_rmtlr` predicts immune response using regularized multi-task
#' linear algorithm.
#' This algorithm employs model parameters learned during training on different
#' types of data in order to compute the immune response.
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param view_name input view name
#' @param view_info input view information of its composition.
#' @param view_data input view data as a list. Each item of the list corresponds
#' to a certain view.
#' @param learned_model parameters learned during training with cross-validation.
#' @param verbose A logical value indicating whether to display messages about the prediction process.
#'
#' @return A matrix with the predictions obtained by applying the model on the
#' view input data
#'
#' @examples
#' # TODOTODO
predict_with_rmtlr <- function(view_name,
                               view_info,
                               view_data,
                               learned_model,
                               verbose = TRUE) {

  # Initialize variables
  P <- length(view_info)
  K <- 100
  standardize_any <- TRUE
  models <- names(learned_model[[1]]$model$cv.glmnet.features)
  drugs <- colnames(learned_model[[1]]$model$cv.glmnet.features[[1]])
  Ndrug <- length(drugs)
  predictions <- predictions_all_tasks <- predictions_all_models <- list()

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

  # Per task, per view
  predictions[[view_name]] <- matrix(NA,
    nrow = nrow(view_data_new[[1]]), ncol = K,
    dimnames = list(rownames(view_data_new[[1]]), seq(1, K, 1))
  )

  predictions_all_models <- do.call(c, lapply(models, function(X) {
    predictions_all_models[[X]] <- predictions
    return(predictions_all_models)
  }))

  predictions_all_tasks <- do.call(c, lapply(drugs, function(X) {
    predictions_all_tasks[[X]] <- predictions_all_models
    return(predictions_all_tasks)
  }))

  for (i in 1:K) {
    state <- learned_model[[i]]$model$cv.glmnet.features
    features_learning <- lapply(1:length(view_info), function(x) {
      names(learned_model[[i]]$mas.mea.learning.X[[x]])
    })
    prediction_X <- view_data_new

    # Display progress bar:
    if (verbose){
      width <- options()$width
      cat(paste0(rep("=", i / K * width), collapse = ""))
      Sys.sleep(.05)
      if (i == K) {
        cat("\n")
      } else {
        cat(" \r")
      }
    }

    # standardize
    if (standardize_any == TRUE) {
      for (m in 1:P) {

        # Check features availability
        keep_pos <- stats::na.omit(match(colnames(prediction_X[[m]]), features_learning[[m]]))
        keep_names <- intersect(colnames(prediction_X[[m]]), features_learning[[m]])
        prediction_X[[m]] <- prediction_X[[m]][, keep_names]

        # Normalization should be done taking into account the train set
        learned_model[[i]]$mas.mea.learning.X[[m]] <- learned_model[[i]]$mas.mea.learning.X[[m]][keep_pos]
        learned_model[[i]]$mas.std.learning.X[[m]] <- learned_model[[i]]$mas.std.learning.X[[m]][keep_pos]
        names(learned_model[[i]]$mas.std.learning.X[[m]]) <- names(learned_model[[i]]$mas.mea.learning.X[[m]])

        prediction_X[[m]] <- standardization(
          X = prediction_X[[m]], mean = learned_model[[i]]$mas.mea.learning.X[[m]],
          sd = learned_model[[i]]$mas.std.learning.X[[m]]
        )
      }
    }

    # perform prediction
    prediction_cv <- lapply(state, function(X) {
      rmtlr_test(prediction_X, X)
    })

    # save predictions
    for (X in drugs) {
      predictions_all_tasks[[X]][["1se.mse"]][[view_name]][, i] <- prediction_cv$`1se.mse`[, X]
      predictions_all_tasks[[X]][["min.mse"]][[view_name]][, i] <- prediction_cv$min.mse[, X]
    }
  }
  summary_pred <- predictions_all_tasks
  return(summary_pred)
}
