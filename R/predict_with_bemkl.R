#' Make predictions using BEMKL
#'
#' `predict_with_bemkl` predicts immune response using bayesian efficient
#' multi-kernel algorithm. This algorithm employs model parameters learned
#' during training on different types of data in order to compute the immune response.
#'
#' @importFrom pdist pdist
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
predict_with_bemkl <- function(view_name,
                               view_info,
                               view_data,
                               learned_model,
                               verbose = TRUE) {

  # Initialize variables
  P <- length(view_info)
  K <- 100
  standardize_any <- TRUE
  drugs <- names(learned_model[[1]]$performances$MSE)
  Ndrug <- length(drugs)
  predictions <- predictions_all_tasks <- list()

  # Per task, per view
  predictions[[view_name]] <- matrix(NA,
    nrow = nrow(view_data[[1]]), ncol = K,
    dimnames = list(rownames(view_data[[1]]), seq(1, K, 1))
  )

  predictions_all_tasks <- do.call(c, lapply(drugs, function(X) {
    predictions_all_tasks[[X]] <- predictions
    return(predictions_all_tasks)
  }))

  # Per iteration
  for (i in 1:K) {
    state <- learned_model[[i]]$model
    learning_X <- learned_model[[i]]$training_set
    prediction_X <- lapply(names(view_info), function(x) {
      view_data[[x]]
    })

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
        if (view_info[m] != "jaccard") {
          # Check same features availability
          keep_pos <- na.omit(match(colnames(prediction_X[[m]]), colnames(learning_X[[m]])))
          keep_names <- intersect(colnames(prediction_X[[m]]), colnames(learning_X[[m]]))
          prediction_X[[m]] <- prediction_X[[m]][, keep_names]
          learning_X[[m]] <- learning_X[[m]][, keep_names]

          # Normalization should be done taking into account the train set. #
          learned_model[[i]]$mas.mea.learning.X[[m]] <- learned_model[[i]]$mas.mea.learning.X[[m]][keep_pos]
          learned_model[[i]]$mas.std.learning.X[[m]] <- learned_model[[i]]$mas.std.learning.X[[m]][keep_pos]

          prediction_X[[m]] <- standardization(
            prediction_X[[m]], learned_model[[i]]$mas.mea.learning.X[[m]],
            learned_model[[i]]$mas.std.learning.X[[m]]
          )
        }
      }
    }

    # compute prediction kernel
    Nlearning <- nrow(learning_X[[1]])
    Nprediction <- nrow(prediction_X[[1]])
    Kx_prediction <- array(rep(0, Nlearning * Nprediction * P), c(Nlearning, Nprediction, P))

    for (m in 1:P) {
      Kx_prediction[, , m] <- exp(-(as.matrix(pdist(learning_X[[m]], prediction_X[[m]])))^2 / ncol(learning_X[[m]]) / 2)
    }

    Ktest <- Kx_prediction # should be an Ntra x Ntest x P matrix containing similarity values between training and test samples

    # perform prediction
    prediction <- bemkl_supervised_multioutput_regression_variational_test(Ktest, state)
    predictions <- t(prediction$Y$mu)
    colnames(predictions) <- drugs
    rownames(predictions) <- rownames(prediction_X[[1]])

    # save predictions
    for (X in drugs) {
      predictions_all_tasks[[X]][[view_name]][, i] <- predictions[, X]
    }
  }

  summary_pred <- predictions_all_tasks
  return(summary_pred)
}
