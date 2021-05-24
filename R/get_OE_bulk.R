#' Calculate overall expression (OE)
#'
#' `get_OE_bulk` obtained from literature to calculate Immune resistance program
#' (Jerby-Arnon et al., 2018)
#'
#' @importFrom arules discretize
#'
#' @param r list TODOTODO - needs some more info?
#' @param gene_sign string
#' @param num_rounds integer
#' @param full_flag boolean
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return Random score
#'
#' @examples
#' # TODOTODO
get_OE_bulk <- function(r,
                        gene_sign = NULL,
                        num_rounds = 1000,
                        full_flag = FALSE,
                        verbose = TRUE) {
  set.seed(1234)
  r$genes_mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes_mean, FUN = "-")
  r$genes_dist <- r$genes_mean
  r$genes_dist_q <- arules::discretize(r$genes_dist, n.cat = 50)
  r$sig_scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene_sign))
  sig_names <- names(gene_sign)
  colnames(r$sig_scores) <- sig_names
  r$sig_scores_raw <- r$sig_scores
  rand_flag <- is.null(r$rand_scores) | !all(is.element(names(gene_sign), colnames(r$rand_scores)))
  if (rand_flag) {
    #if (verbose) message("Computing also random scores...", "\n")
    r$rand_scores <- r$sig_scores
  }
  for (i in sig_names) {
    b_sign <- is.element(r$genes, gene_sign[[i]])
    if (sum(b_sign) < 2) {
      next()
    }
    if (rand_flag) {
      rand_scores <- get_semi_random_OE(r, r$genes_dist_q, b_sign, num_rounds = num_rounds)
    } else {
      rand_scores <- r$rand_scores[, i]
    }
    raw_scores <- colMeans(r$zscores[b_sign, ])
    final_scores <- raw_scores - rand_scores
    r$sig_scores[, i] <- final_scores
    r$sig_scores_raw[, i] <- raw_scores
    r$rand_scores[, i] <- rand_scores
  }
  if (full_flag) {
    return(r)
  }
  sig_scores <- r$sig_scores
  return(sig_scores)
}
