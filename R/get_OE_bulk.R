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
#'
#' @return Random score
#'
#' @examples
#' # TODOTODO
get_OE_bulk <- function(r,
                        gene_sign = NULL,
                        num_rounds = 1000,
                        full_flag = FALSE) {
  set.seed(1234)
  r$genes.mean <- rowMeans(r$tpm)
  r$zscores <- sweep(r$tpm, 1, r$genes.mean, FUN = "-")
  r$genes.dist <- r$genes.mean
  r$genes.dist.q <- arules::discretize(r$genes.dist, n.cat = 50)
  r$sig.scores <- matrix(data = 0, nrow = ncol(r$tpm), ncol = length(gene_sign))
  sig.names <- names(gene_sign)
  colnames(r$sig.scores) <- sig.names
  r$sig.scores.raw <- r$sig.scores
  rand.flag <- is.null(r$rand.scores) | !all(is.element(names(gene_sign), colnames(r$rand.scores)))
  if (rand.flag) {
    # TODOTODO: maybe use message instead - it is handled in a more gentle way and could be suppressed in practical manners ;)
    # TODOTODO: could apply to other print commands
    message("Computing also random scores...")
    r$rand.scores <- r$sig.scores
  }
  for (i in sig.names) {
    b.sign <- is.element(r$genes, gene_sign[[i]])
    if (sum(b.sign) < 2) {
      next()
    }
    if (rand.flag) {
      rand.scores <- get_semi_random_OE(r, r$genes.dist.q, b.sign, num_rounds = num_rounds)
    } else {
      rand.scores <- r$rand.scores[, i]
    }
    raw.scores <- colMeans(r$zscores[b.sign, ])
    final.scores <- raw.scores - rand.scores
    r$sig.scores[, i] <- final.scores
    r$sig.scores.raw[, i] <- raw.scores
    r$rand.scores[, i] <- rand.scores
  }
  if (full_flag) {
    return(r)
  }
  sig.scores <- r$sig.scores
  return(sig.scores)
}
