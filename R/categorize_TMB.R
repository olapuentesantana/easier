#' Convert tumor mutational burden (TMB) to categorical.
#'
#' `categorize_TMB` encodes tumor mutational variable (TMB) numerical variable to categorical.
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @param TMB A numeric vector with tumor mutational burden values.
#' @param thresholds A numeric vector to specify thresholds to be used. Default thresholds are low (<100), moderate (100-400) and high TMB (>400).
#'
#' @return A numeric vector assigning each sample a class from 1 to 3.
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#'data(cds)
#'mariathasan_data <- preprocess_mariathasan(cds)
#'rm(cds)
#'
#'# retrieve Tumor Mutational Burden
#'TumorMutationalBurden <- clinical_data[, "FMOne mutation burden per MB"]
#'names(TumorMutationalBurden) <- rownames(clinical_data)
#'
#' # Convert TMB continous values into categories
#' TMB_cat <- categorize_TMB(TMB = TumorMutationalBurden)
categorize_TMB <- function(TMB, thresholds=NULL) {
  vTert <- stats::quantile(TMB , c(0:3/3))
  if (is.null(thresholds)){
    x_out <- cut(TMB, vTert, include.lowest = TRUE, labels = c(1, 2, 3))
  }else if ((is.numeric(thresholds)) & (length(thresholds)==2)){
    # x_out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)), include.lowest = T, labels = c(1, 2, 3))
    TMB_out <- cut(TMB, c(min(TMB, na.rm = TRUE), 100, 400, max(TMB, na.rm = TRUE)), include.lowest = FALSE, labels = c(1, 2, 3))
  }
  TMB_out <- as.numeric(TMB_out)
  names(TMB_out) <- names(TMB)
  return(TMB_out)
}
