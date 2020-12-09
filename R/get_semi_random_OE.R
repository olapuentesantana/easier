#' Calculate random scores
#'
#' \code{get_semi_random_OE} obtained from literature to calculate Immune resistance program  (Jerby-Arnon et al., 2018)
#'
#' @param r list
#' @param genes.dist.q integer
#' @param b.sign boolean
#' @param num.rounds integer
#' @param full.flag boolean
#'
#' @return Random score
#'
#' @€xamples
#' # TODOTODO
get_semi_random_OE <- function(r,genes.dist.q,b.sign,num.rounds = 1000,full.flag = F){
  # Previous name: get.random.sig.scores
  sign.q <- as.matrix(table(genes.dist.q[b.sign]))
  q <- rownames(sign.q)
  idx.all <- c()
  B <- matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q <- matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
  for (i in 1:nrow(sign.q)){
    num.genes <- sign.q[i]
    if(num.genes > 0){
      idx <- which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj <- sample(idx,num.genes)
        Q[i,j] <- sum(B[idxj,j]==T)
        B[idxj,j] <- T
      }
    }
  }
  rand.scores <- apply(B,2,function(x) colMeans(r$zscores[x,]))
  if(full.flag){return(rand.scores)}
  rand.scores <- rowMeans(rand.scores)
  return(rand.scores)
}
