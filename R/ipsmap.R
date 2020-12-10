#' Calculate Immunophenoscore value
#'
#' \code{ipsmap} obtained from literature to calculate immunophenoscore.  (Charoentong et al., 2017)
#'
#' @param x numeric value per sample, needed to check more precisely.
#'
#' @return IPS per sample
#'
#' @examples
#' # TODOTODO
ipsmap <- function (x) {
  if (x <= 0) {
    ips <- 0
  } else {
    if (x >= 3) {
      ips <- 10
    } else {
      ips <- round(x * 10/ 3, digits = 0)
    }
  }
  return(ips)
}
