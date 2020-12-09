#' Compute cytolytic activity score
#'
#' \code{compute_CYT} computes cytolytic activity score as the geometric mean of immune cytolytic genes
#' (Rooney et al., 2015).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#'
#' @€xamples
#' # TODOTODO
compute.CYT <- function(RNA.tpm){

  # Literature genes
  CYT.read <- c("GZMA", "PRF1")
  match_CYT.genes <- match(CYT.read, rownames(RNA.tpm))

  if (anyNA(match_CYT.genes)){
    warning(paste0("differenty named or missing signature genes : \n", paste(CYT.read[!CYT.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_CYT.genes <- stats::na.omit(match_CYT.genes)
  }

  # Subset RNA.tpm
  subset_RNA.tpm <- RNA.tpm[match_CYT.genes, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}
