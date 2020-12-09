#' Compute T cell-inflamed signature score
#'
#' \code{compute_ayersTcellInfl} computes T cell-inflamed signature score by taking a weighted sum of
#'  the housekeeping normalized values of the T cell-inflamed signature genes
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=T cell-inflamed signature score
#'
#' @export
#'
#' @€xamples
#' # TODOTODO
compute.Tcell_inflamed <- function(RNA.tpm){

  # Literature genes
  Tcell_inflamed.read <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                                  "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT")
  Housekeeping.read <- c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2", "UBB", "TBP", "SDHA") # C14orf102 = NRDE2
  weights <- data.frame(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                        CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                        PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767, check.names = FALSE)

  # Some genes might have other name: case for "C14orf102", it's called "NRDE2", be careful
  if (any(rownames(RNA.tpm) %in% "C14orf102")){
    cat("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "C14orf102"] <- "NRDE2"
  }

  match_genes.housekeeping <- match(Housekeeping.read, rownames(RNA.tpm))
  match_genes.predictors <- match(Tcell_inflamed.read, rownames(RNA.tpm))

  if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
    tmp <- c(Tcell_inflamed.read, Housekeeping.read)
    warning(c("differenty named or missing signature genes : \n", paste(tmp[!tmp %in% rownames(RNA.tpm)], collapse = "\n")))
    match_genes.housekeeping <- stats::na.omit(match_genes.housekeeping)
    match_genes.predictors <- stats::na.omit(match_genes.predictors)
  }

  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  ## housekeeping
  log2.RNA.tpm.housekeeping <- log2.RNA.tpm[match_genes.housekeeping, ]
  ## predictors
  log2.RNA.tpm.predictors <- log2.RNA.tpm[match_genes.predictors, ]
  weights <- weights[,rownames(log2.RNA.tpm.predictors)]

  # Housekeeping normalization
  average.log2.RNA.tpm.housekeeping <- apply(log2.RNA.tpm.housekeeping, 2, mean)
  log2.RNA.tpm.predictors.norm <- sweep(log2.RNA.tpm.predictors, 2, average.log2.RNA.tpm.housekeeping, FUN = "-")

  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.RNA.tpm.predictors.norm), colnames(as.vector(weights)))
  score <- t(log2.RNA.tpm.predictors.norm[tidy,]) %*% t(as.vector(weights))

  message("Tcell_inflamed score computed")
  return(data.frame( Tcell_inflamed = score, check.names = FALSE))
}
