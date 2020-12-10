#' Multi-task elastic net algorithm implements the model training coefficients and
#' hyperparameters on the test set.
#'
#' \code{multi_task_EN_test} implements mutli-task elastic net regression on the test data
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param x.test A matrix with data required to perform predictions.
#' @param coef.matrix A coefficient matrix obtained during training the model (rows = features and columns = tasks).
#'
#' @return A matrix with predictions for each sample and task.
#'
#' @examples
#' # TODOTODO
multi_task_EN_test <- function(x.test, coef.matrix){

  # Keep intercept
  Intercept <- as.matrix(coef.matrix)[1,]

  # Remove intercept from coef matrix
  coef <- as.matrix(coef.matrix)[-1, ,drop = FALSE]

  # Combine views
  x.test.combo <- do.call(cbind, lapply(1:length(x.test), function(x){tmp = x.test[[x]]}))

  # match features properly
  pos <- stats::na.omit(match(colnames(x.test.combo), rownames(coef)))
  coef <- coef[pos,, drop = FALSE]

  if (length(coef) > 1){
    Slope <-  coef
    rownames(Slope) <- gsub(" ","",rownames(Slope)) # In case of Dorothea is needed, do not affect the other features
    fit.pred <- t(matrix(as.matrix(Intercept), nrow = ncol(Slope), ncol = nrow(x.test.combo))
                  + t(Slope) %*% t(as.matrix(x.test.combo[, rownames(Slope)])))
  }else{
    Slope <- 0
    fit.pred <- matrix(as.matrix(Intercept), nrow = ncol(Slope), ncol = nrow(x.test[[1]]))
  }
  return(fit.pred)
}
