#' Regularized multi-task linear regression algorithm implements the model
#' training coefficients and hyperparameters on the test set.
#'
#' `rmtlr_test` implements regularized multi-task linear regression on the test data
#'
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param x_test A matrix with data required to perform predictions.
#' @param coef_matrix A coefficient matrix obtained during training the model
#' (rows = features and columns = tasks).
#'
#' @return A matrix with predictions for each sample and task.
#'
#' @examples
#' # TODOTODO
rmtlr_test <- function(x_test,
                       coef_matrix) {

  # Keep intercept
  intercept <- as.matrix(coef_matrix)[1, ]

  # Remove intercept from coef matrix
  coef <- as.matrix(coef_matrix)[-1, , drop = FALSE]

  # Combine views
  x_test_combo <- do.call(cbind, lapply(1:length(x_test), function(x) {
    tmp <- x_test[[x]]
  }))

  # match features properly
  pos <- stats::na.omit(match(colnames(x_test_combo), rownames(coef)))
  coef <- coef[pos, , drop = FALSE]

  if (length(coef) > 1) {
    slope <- coef
    fit_pred <- t(matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test_combo))
    + t(slope) %*% t(as.matrix(x_test_combo[, rownames(slope)])))
  } else {
    slope <- 0
    fit_pred <- matrix(as.matrix(intercept), nrow = ncol(slope), ncol = nrow(x_test[[1]]))
  }
  return(fit_pred)
}
