#' Compute immunophenoscore
#'
#' `compute_IPS` computes the immunophenoscore using source code provided by
#' original publication (Charoentong et al., 2017).
#'
#' @importFrom stats sd
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=Immunophenoscore
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' IPS <- compute_IPS(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(IPS)
compute_IPS <- function(RNA_tpm,
                        verbose = TRUE) {

  # Log2 transformation:
  log2_RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))
  sample_names <- colnames(log2_RNA_tpm)

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
  GVEC <- row.names(log2_RNA_tpm)

  # Genes names in IPS genes file
  VEC <- as.vector(IPSG$GENE)

  # Match IPS genes with genes in expression file
  ind <- which(is.na(match(VEC, GVEC)))

  # List genes missing or differently named
  MISSING_GENES <- VEC[ind]
  dat <- IPSG[ind, ]
  if (length(MISSING_GENES) > 0) {
    warning("differently named or missing genes : \n", MISSING_GENES, "\n")
    for (x in 1:length(ind)) {
      print(IPSG[ind, ])
    }
  }

  # calculation
  for (i in 1:length(sample_names)) {
    GE <- log2_RNA_tpm[[i]]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (log2_RNA_tpm[as.vector(IPSG$GENE), i] - mGE) / sGE
    W1 <- IPSG$WEIGHT
    WEIGHT <- NULL
    MIG <- NULL
    k <- 1
    for (gen in unique_ips_genes) {
      MIG[k] <- mean(Z1[which(as.vector(IPSG$NAME) == gen)], na.rm = TRUE)
      WEIGHT[k] <- mean(W1[which(as.vector(IPSG$NAME) == gen)])
      k <- k + 1
    }
    WG <- MIG * WEIGHT
    MHC[i] <- mean(WG[1:10])
    CP[i] <- mean(WG[11:20])
    EC[i] <- mean(WG[21:24])
    SC[i] <- mean(WG[25:26])
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    score[i] <- ipsmap(AZ[i])
  }
  names(score) <- sample_names

  if (verbose) message("IPS score computed")
  return(data.frame(IPS = score, check.names = FALSE))
}
