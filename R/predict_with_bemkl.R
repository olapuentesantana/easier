#' Make predictions using BEMKL
#'
#' \code{predict_with_bemkl} predicts immune response using bayesian efficient multi-kernel algorithm.
#' This algorithm employs model parameters learned during training on different types of data in
#' order to compute the immune response.
#'
#' @importFrom pdist pdist
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param view_name input view name
#' @param view_info input view information of its composition.
#' @param view_data input view data as a list. Each item of the list correponds to a certain view.
#' @param learned_model parameters learned during training with cross-validation.
#'
#' @return A matrix with the predictions obtained by applying the model on the view input data
#'
#' @examples
#' # TODOTODO
predict_with_bemkl <- function(view_name, view_info, view_data, learned_model){

  # Initialize variables
  P <- length(view_info)
  K <- 100
  standardize_any <- TRUE
  drugs <- names(learned_model[[1]]$performances$MSE)
  Ndrug <- length(drugs)
  predictions <- predictions.all.tasks <- list()

  # Per task, per view
  predictions[[view_name]] <- matrix(NA, nrow = nrow(view_data[[1]]), ncol = K,
                                     dimnames = list(rownames(view_data[[1]]), seq(1, K, 1)))

  predictions.all.tasks <- do.call(c, lapply(drugs, function(X){
    predictions.all.tasks[[X]] <- predictions
    return(predictions.all.tasks)
  }))

  # Per iteration
  for (i in 1:K){
    state <- learned_model[[i]]$model
    learning.X <- learned_model[[i]]$training_set
    prediction.X <- lapply(names(view_info), function(x){view_data[[x]]})

    message("Iteration ", i,"\n")
    # standardize
    if (standardize_any==TRUE){
      for (m in 1:P){

        if (view_info[m] != "jaccard"){
          # Check same features availability
          keep_pos <- na.omit(match(colnames(prediction.X[[m]]), colnames(learning.X[[m]])))
          keep_names <- intersect(colnames(prediction.X[[m]]), colnames(learning.X[[m]]))
          prediction.X[[m]] <- prediction.X[[m]][,keep_names]
          learning.X[[m]] <- learning.X[[m]][,keep_names]

          # Normalization should be done taking into account the train set. #
          learned_model[[i]]$mas.mea.learning.X[[m]] <- learned_model[[i]]$mas.mea.learning.X[[m]][keep_pos]
          learned_model[[i]]$mas.std.learning.X[[m]] <- learned_model[[i]]$mas.std.learning.X[[m]][keep_pos]

          prediction.X[[m]] <- standarization(prediction.X[[m]], learned_model[[i]]$mas.mea.learning.X[[m]],
                                              learned_model[[i]]$mas.std.learning.X[[m]])

        }
      }
    }

    # compute prediction kernel
    Nlearning <- nrow(learning.X[[1]])
    Nprediction <- nrow(prediction.X[[1]])
    Kx_prediction <- array(rep(0, Nlearning*Nprediction*P), c(Nlearning, Nprediction, P))

    for (m in 1:P){
      Kx_prediction[, , m] <- exp(-(as.matrix(pdist(learning.X[[m]], prediction.X[[m]]))) ^ 2 / ncol(learning.X[[m]]) / 2)
    }

    Ktest <- Kx_prediction #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples

    #perform prediction
    prediction <- bemkl_supervised_multioutput_regression_variational_test(Ktest, state)
    predictions <- t(prediction$Y$mu)
    colnames(predictions) <- drugs
    rownames(predictions) <- rownames(prediction.X[[1]])

    # save predictions
    for (X in drugs){
      predictions.all.tasks[[X]][[view_name]][,i] <- predictions[,X]
    }
  }

  summary_pred <- predictions.all.tasks
  return(summary_pred)
}
