#' Make predictions using multi-task elastic net
#'
#' `predict_with_multitaskelasticnet` predicts immune response using multi-task
#' elastic net algorithm.
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
#'
#' @return A matrix with the predictions obtained by applying the model on the
#' view input data
#'
#' @examples
#' # TODOTODO
predict_with_multitaskelasticnet <- function(view_name,
                                             view_info,
                                             view_data,
                                             learned_model) {

  # Initialize variables
  P <- length(view_info)
  K <- 100
  standardize_any <- TRUE
  models <- names(learned_model[[1]]$model$cv.glmnet.features)
  drugs <- colnames(learned_model[[1]]$model$cv.glmnet.features[[1]])
  Ndrug <- length(drugs)
  predictions <- predictions.all.tasks <- predictions.all.models <- list()

  # Algorithm do not deal with NA values: here we removed features with NA values, patients with all NA values are not removed
  view_data_new <- lapply(names(view_data), function(x){
    tmp_data <- view_data[[x]]
    if (all(is.na(apply(tmp_data, 1, sum)))){
      NA_sum <- apply(tmp_data, 2, sum)
      tmp_data <- tmp_data[,!is.na(NA_sum)]
    }
    return(tmp_data)
  })
  names(view_data_new) <- names(view_data)

  # Per task, per view
  predictions[[view_name]] <- matrix(NA, nrow = nrow(view_data_new[[1]]), ncol = K,
                                     dimnames = list(rownames(view_data_new[[1]]), seq(1, K, 1)))

  predictions.all.models <- do.call(c, lapply(models, function(X){
    predictions.all.models[[X]] <- predictions
    return(predictions.all.models)
  }))

  predictions.all.tasks <- do.call(c, lapply(drugs, function(X){
    predictions.all.tasks[[X]] <- predictions.all.models
    return(predictions.all.tasks)
  }))

  for (i in 1:K){

    state <- learned_model[[i]]$model$cv.glmnet.features
    features.learning <- lapply(1:length(view_info), function(x){names(learned_model[[i]]$mas.mea.learning.X[[x]])})
    prediction.X <- view_data_new

    # Display progress bar:
    width <- options()$width
    cat(paste0(rep('=', i / K * width), collapse = ''))
    Sys.sleep(.05)
    if (i == K) cat('\n')
    else cat(' \r')

    # standardize
    if (standardize_any==TRUE){
      for (m in 1:P){

        # Check features availability
        keep_pos <- stats::na.omit(match(colnames(prediction.X[[m]]), features.learning[[m]]))
        keep_names <- intersect(colnames(prediction.X[[m]]), features.learning[[m]])
        prediction.X[[m]] <- prediction.X[[m]][,keep_names]

        # Normalization should be done taking into account the train set
        learned_model[[i]]$mas.mea.learning.X[[m]] <- learned_model[[i]]$mas.mea.learning.X[[m]][keep_pos]
        learned_model[[i]]$mas.std.learning.X[[m]] <- learned_model[[i]]$mas.std.learning.X[[m]][keep_pos]
        names(learned_model[[i]]$mas.std.learning.X[[m]]) <- names(learned_model[[i]]$mas.mea.learning.X[[m]])

        prediction.X[[m]] <- standardization(X = prediction.X[[m]], mean = learned_model[[i]]$mas.mea.learning.X[[m]],
                                            sd = learned_model[[i]]$mas.std.learning.X[[m]])
      }
    }

    # perform prediction
    prediction_cv <- lapply(state, function(X){multi_task_EN_test(prediction.X, X)})

    # save predictions
    for (X in drugs){
      predictions.all.tasks[[X]][["1se.mse"]][[view_name]][,i] <- prediction_cv$`1se.mse`[,X]
      predictions.all.tasks[[X]][["min.mse"]][[view_name]][,i] <- prediction_cv$min.mse[,X]
    }
  }
  summary_pred <- predictions.all.tasks
  return(summary_pred)
}
