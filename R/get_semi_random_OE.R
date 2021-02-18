#' Calculate random scores
#'
#' `get_semi_random_OE` obtained from literature to calculate Immune resistance
#' program  (Jerby-Arnon et al., 2018) TODOTODO - out of sync?
#'
#' @param r list TODOTODO - more info here as well
#' @param genes_dist_q integer
#' @param b_sign boolean
#' @param num_rounds integer
#' @param full_flag boolean
#'
#' @return Random score
#'
#' @examples
#' # TODOTODO
get_semi_random_OE <- function(r,
                               genes_dist_q,
                               b_sign,
                               num_rounds = 1000,
                               full_flag = FALSE) {
  # Previous name: get.random.sig.scores
  sign.q <- as.matrix(table(genes_dist_q[b_sign]))
  q <- rownames(sign.q)
  idx.all <- c()
  B <- matrix(data = FALSE,nrow = length(genes_dist_q),ncol = num_rounds)
  Q <- matrix(data = 0,nrow = length(genes_dist_q),ncol = num_rounds)
  for (i in 1:nrow(sign.q)){
    num.genes <- sign.q[i]
    if(num.genes > 0){
      idx <- which(is.element(genes_dist_q,q[i]))
      for (j in 1:num_rounds){
        idxj <- sample(idx,num.genes)
        Q[i,j] <- sum(B[idxj,j]==TRUE)
        B[idxj,j] <- TRUE
      }
    }
  }
  rand.scores <- apply(B,2,function(x) colMeans(r$zscores[x,]))
  if(full_flag){return(rand.scores)}
  rand.scores <- rowMeans(rand.scores)
  return(rand.scores)
}
