#' Define tumor mutational burden (TMB) categories.
#'
#' Encodes tumor mutational burden (TMB)
#' from numerical variable into categorical.
#'
#' @importFrom stats quantile
#'
#' @param TMB numeric vector with tumor mutational burden values.
#' @param thresholds numeric vector to specify thresholds to be used.
#' Default thresholds are low (<100), moderate (100-400) and
#' high TMB (>400).
#'
#' @return A numeric vector assigning each sample a class
#' from 1 (low TMB) to 3 (high TMB).
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
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#'
#' TMB <- colData(dataset_mariathasan)[["TMB"]]
#' names(TMB) <- colData(dataset_mariathasan)[["pat_id"]]
#'
#' # Convert TMB continous values into categories
#' TMB_cat <- categorize_TMB(TMB = TMB)
categorize_TMB <- function(TMB, thresholds = NULL) {
    vTert <- stats::quantile(TMB, probs = c(0:3 / 3), na.rm = TRUE)
    if (is.null(thresholds)) {
        TMB_out <- cut(TMB, vTert, include.lowest = TRUE, labels = c(1, 2, 3))
    } else if ((is.numeric(thresholds)) & (length(thresholds) == 2)) {
        # x_out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)),
        #             include.lowest = T, labels = c(1, 2, 3))
        TMB_out <- cut(TMB, c(
            min(TMB, na.rm = TRUE), 100, 400,
            max(TMB, na.rm = TRUE)
        ),
        include.lowest = FALSE, labels = c(1, 2, 3)
        )
    }
    TMB_out <- as.numeric(TMB_out)
    names(TMB_out) <- names(TMB)
    return(TMB_out)
}
