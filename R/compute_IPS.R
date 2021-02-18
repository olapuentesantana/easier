#' Compute immunophenoscore
#'
#' `compute_IPS` computes the immunophenoscore using source code provided by
#' original publication (Charoentong et al., 2017).
#'
#' @importFrom stats sd
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Immunophenoscore
#'
#' @export
#'
#' @examples
#' # TODOTODO
compute.IPS <- function(RNA.tpm) {

  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))
  sample_names <- colnames(log2.RNA.tpm)

  # Literature genes and corresponding weights
  IPSG <- IPSG_read
  unique_ips_genes <- as.vector(unique(IPSG$NAME))

  score <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL

  # Gene names in expression file
  GVEC <- row.names(log2.RNA.tpm)

  # Genes names in IPS genes file
  VEC <- as.vector(IPSG$GENE)

  # Match IPS genes with genes in expression file
  ind <- which(is.na(match(VEC, GVEC)))

  # List genes missing or differently named
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind, ]
  if (length(MISSING_GENES) > 0) {
    cat("differently named or missing genes: ", MISSING_GENES,"\n")
    for (x in 1:length(ind)) {
      print(IPSG[ind, ])
    }
  }

  # calculation
  for (i in 1:length(sample_names)) {
    GE <- log2.RNA.tpm[[i]]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (log2.RNA.tpm[as.vector(IPSG$GENE),i] - mGE) / sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    for (gen in unique_ips_genes) {
      MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)], na.rm=TRUE)
      WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k <- k+1
    }
    WG <- MIG*WEIGHT
    MHC[i] <- mean(WG[1:10])
    CP[i] <- mean(WG[11:20])
    EC[i] <- mean(WG[21:24])
    SC[i] <- mean(WG[25:26])
    AZ[i] <- sum(MHC[i],CP[i],EC[i],SC[i])
    score[i] <- ipsmap(AZ[i])
  }
  names(score) <- sample_names

  message("IPS score computed")
  return(data.frame(IPS = score, check.names = FALSE))
}
