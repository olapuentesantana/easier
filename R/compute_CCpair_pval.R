#' compute_CCpair_pval
#'
#' \code{compute_CCpair_pval} compute corresponding p-values by comparing the score with the null hypothesis computed for each cell-cell pair.
#'
#' @export
#'
#' @param value CC pair score
#' @param null.hypothesis hypothesized distribution for each cell-cell pair
#'
#' @return p-value for each cell-cell pair
#'
#' @examples
#' # TODOTODO
compute_CCpair_pval <- function(value, null.hypothesis){

  sum(null.hypothesis >= value)/length(null.hypothesis)

}
